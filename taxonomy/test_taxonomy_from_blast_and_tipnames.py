from taxonomy_from_blast_and_tipnames import (coreset_id_to_rdp_taxa,
        each_query_top_coresets, get_query_common_lineages,
        CORESET_RDP_NDS_PATH)
from py_util.unit_test import TestCase, main, set_trace, StringIO

class Test_coreset_id_to_rdp_taxa(TestCase):
    def test_basic(self):
        res = coreset_id_to_rdp_taxa(StringIO(coreset_rdp_nds_test))
        self.assertEqual(res, coreset_rdp_exp)

class Test_each_query_top_coresets(TestCase):
    def test_basic(self):
        res = each_query_top_coresets(StringIO(blast_test))
        self.assertEqual(list(res), query_top_coresets_exp)

class Test_get_query_common_lineages(TestCase):
    def test_basic(self):
        res = get_query_common_lineages(StringIO(blast_test),
                file(CORESET_RDP_NDS_PATH))
        self.assertEqual(res, common_lineage_exp)


### test data
coreset_rdp_nds_test = """\
63788	Bacteria; Genera_incertae_sedis_WS3; WS3
63696	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Pseudonocardineae; Actinosynnemataceae; Saccharothrix
63834	Bacteria; Firmicutes; Bacilli; Bacillales; Bacillaceae; unclassified_Bacillaceae
00000	   
"""

coreset_rdp_exp = {
'63696': ['bacteria', 'actinobacteria', 'actinobacteria', 'actinomycetales',
        'actinosynnemataceae'], #subclass ignored
 '63788': ['bacteria', 'genera_incertae_sedis_ws3'], #unknown class ignored
 '63834': ['bacteria', 'firmicutes', 'bacilli', 'bacillales', 'bacillaceae'],
 # empty taxon ignored
 }

blast_test = """\
# BLASTN 2.2.15 [Oct-15-2006]
# Query: C20_b12_clip_784_926 C20_c20_clip_784_926
# Database: greengenes_coreset_x
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
C20_b12_clip_784_926	gi|67942541|greengenes|131682	98.61	144	2	0	1	144	755	898	9e-73	 270
C20_b12_clip_784_926	gi|67942756|greengenes|133374	98.61	144	2	0	1	144	762	905	9e-73	 270
C20_b12_clip_784_926	gi|71063920|greengenes|132359	96.53	144	5	0	1	144	669	812	1e-65	 246

# BLASTN 2.2.15 [Oct-15-2006]
# Query: B936_clip_784_926
# Database: greengenes_coreset_x
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
B936_clip_784_926	gi|13537348|greengenes|52643	99.29	140	1	0	1	140	771	910	9e-73	 270
B936_clip_784_926	gi|18693131|greengenes|25165	99.26	136	1	0	1	136	774	909	2e-70	 262
"""
query_top_coresets_exp = [
('C20_b12_clip_784_926', ['131682', '133374']), #query is the first word
('B936_clip_784_926', ['52643'])]

common_lineage_exp = { #not checked carefully
 'B936_clip_784_926': ['bacteria',
                       'bacteroidetes',
                       'bacteroidetes',
                       'bacteroidales',
                       'bacteroidaceae'],
 'C20_b12_clip_784_926': ['bacteria',
                          'firmicutes',
                          'clostridia',
                          'clostridiales']}
if __name__ == '__main__':
    main()
