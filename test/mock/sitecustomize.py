import ihm
import os

# If we're running from an SGE job, override the from_pubmed_id() function
# to return a cached value, since we don't have network access (needed to
# query PubMed directly)

def mock_from_pubmed(cls, pubmed_id):
    return ihm.Citation(
            pmid=29539637,
            title='Integrative structure and functional anatomy of a '
                  'nuclear pore complex.',
            volume=555, page_range=(475,482), year=2018, authors=[
                'Kim SJ', 'Fernandez-Martinez J', 'Nudelman I', 'Shi Y',
                'Zhang W', 'Raveh B', 'Herricks T', 'Slaughter BD', 'Hogan JA',
                'Upla P', 'Chemmama IE', 'Pellarin R', 'Echeverria I',
                'Shivaraju M', 'Chaudhury AS', 'Wang J', 'Williams R',
                'Unruh JR', 'Greenberg CH', 'Jacobs EY', 'Yu Z',
                'de la Cruz MJ', 'Mironska R', 'Stokes DL', 'Aitchison JD',
                'Jarrold MF', 'Gerton JL', 'Ludtke SJ', 'Akey CW', 'Chait BT',
                'Sali A', 'Rout MP'], doi='10.1038/nature26003')

if 'JOB_ID' in os.environ:
    ihm.Citation.from_pubmed_id = classmethod(mock_from_pubmed)
