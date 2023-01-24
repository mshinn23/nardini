__all__ = [
    'DEFAULT_RANDOM_SEED',
    'DEFAULT_PREFIX_NAME',
    'NUM_SCRAMBLED_SEQUENCES',
    'MAPPING_8x8',
    'MAPPING_9x9',
    'TYPEALL_8x8',
    'TYPEALL_9x9',
    'TYPEALL',
    'LABELS_8x8',
    'LABELS_9x9'
]

DEFAULT_RANDOM_SEED     = 0
DEFAULT_PREFIX_NAME     = 'fasta'
NUM_SCRAMBLED_SEQUENCES = 100000


pol  = ['S','T','N','Q','C','H']
pol9 = ['S','T','N','Q','C']
hyd  = ['I','L','M','V']
pos  = ['R','K']
neg  = ['E','D']
aro  = ['F','W','Y']
ala  = ['A']
pro  = ['P']
gly  = ['G']
his  = ['H']


MAPPING_8x8 = {
    'POLAR':        pol,
    'HYDROPHOBIC':  hyd,
    'POSITIVE':     pos,
    'NEGATIVE':     neg,
    'AROMATIC':     aro,
    'ALANINE':      ala,
    'PROLINE':      pro,
    'GLYCINE':      gly,
}

MAPPING_9x9 = {
    'POLAR':        pol9,
    'HYDROPHOBIC':  hyd,
    'POSITIVE':     pos,
    'NEGATIVE':     neg,
    'AROMATIC':     aro,
    'ALANINE':      ala,
    'PROLINE':      pro,
    'GLYCINE':      gly,
    'HISTIDINE':    his,
}

# We use a tuple since this should be immutable
TYPEALL_8x8 = (pol, hyd, pos, neg, aro, ala, pro, gly)
TYPEALL_9x9 = (pol9, hyd, pos, neg, aro, ala, pro, gly, his)
TYPEALL     = TYPEALL_8x8

LABELS_8x8 = ['µ', 'h', '+', '-', 'π', 'A', 'P', 'G']
LABELS_9x9 = ['µ', 'h', '+', '-', 'π', 'A', 'P', 'G', 'H']