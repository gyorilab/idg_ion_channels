import os
import json
import numpy
import boto3
import pickle
import pandas
from indra.statements import stmts_to_json
from indra.assemblers.tsv import TsvAssembler
from indra.assemblers.html import HtmlAssembler
from indra.tools import assemble_corpus as ac
from indra.belief import BeliefEngine
from indra.sources.indra_db_rest import get_statements


def get_channel_statements(channels, ev_limit=100):
    """Get all statements from the database for a list of gene symbols."""
    all_statements = {}
    for channel in channels:
        idbp = get_statements(agents=[channel], ev_limit=ev_limit,
                              best_first=False)
        source_counts = idbp.get_source_counts()
        stmts = filter_out_medscan(idbp.statements, source_counts)
        all_statements[channel] = stmts
    return all_statements


def non_medscan_evidence(stmt, source_counts):
    counts = source_counts.get(stmt.get_hash())
    ev_count = sum(c for k, c in counts.items() if k != 'medscan')
    return ev_count


def filter_out_medscan(stmts, source_counts):
    new_stmts = []
    for stmt in stmts:
        new_evidence = []
        for ev in stmt.evidence:
            if ev.source_api == 'medscan':
                continue
            new_evidence.append(ev)
        if not non_medscan_evidence(stmt, source_counts):
            continue
        new_stmts.append(stmt)
    return new_stmts


def print_statistics(statements):
    counts = sorted([(k, len(s)) for k, s in statements.items()],
                    key=lambda x: x[1], reverse=True)
    raw_counts = [c[1] for c in counts]
    missing = [c[0] for c in counts if c[1] == 0]
    print(f'No statements for channels: {", ".join(missing)}')
    print(f'{counts[0][1]} statements for the top channel {counts[0][0]}')
    print(f'{numpy.mean(raw_counts)} statements on average per channel')


if __name__ == '__main__':
    # Get all channel Statements
    fname = 'data/IDG_target_final.csv'
    df = pandas.read_csv(fname)
    print('Read a total of %d rows from %s' % (len(df), fname))
    df = df[df['idgFamily'] == 'Ion Channel']
    print('Filtered to %d ion channels' % len(df))
    df = df[df['idgTarget'] == True]
    print('Filtered to %d dark ion channels' % len(df))
    channel_gene_names = sorted(list(df['gene']))
    stmts = get_channel_statements(channel_gene_names)
    with open('dark_ion_channel_stmts_v1.pkl', 'wb') as fh:
        pickle.dump(stmts, fh)
    print_statistics(stmts)
