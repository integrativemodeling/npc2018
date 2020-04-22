# Create the accession map used by test/mock/sitecustomize.py

import ihm.reference

codes = [ 'P52891', 'P46673', 'P35729', 'P36161', 'P49687', 'P53011', 'Q04491',
          'Q02647', 'P40368', 'P40477', 'P14907', 'P34077', 'Q02199', 'P48837',
          'P40064', 'P38181', 'P52593', 'P47054', 'Q03790', 'Q05166', 'P32500',
          'Q12445', 'P39685', 'Q02629', 'Q02630', 'P49686', 'Q12315', 'P49687',
          'P20676', 'P39705', 'Q02455', 'P40457' ]

def pp(s):
    indent = 8
    width = 66
    def get_lines(s):
        for i in range(0, len(s), width):
            yield ' ' * indent + "'" + s[i:i+width] + "'"
    return '\n'.join(l for l in get_lines(s))

for code in codes:
    u = ihm.reference.UniProtSequence.from_accession(code)
    print("    '%s': {'db_code':'%s', 'sequence':\n%s},"
          % (code, u.db_code, pp(u.sequence)))
