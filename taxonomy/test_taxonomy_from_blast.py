from taxonomy_from_blast import *
from py_util.unit_test import TestCase, main, StringIO, set_trace
import os


class TaxonomyFromBlast_Tests(TestCase):
    def setUp(self):
        self.tb = TaxonomyFromBlast()
        self.abc = ['a', 'ab', 'abc']
        self.acb =  ['a', 'ac', 'acb']
    def test_init(self): pass
    def test_get_common_lineages(self):
        abc, acb = self.abc, self.acb
        tb = self.tb
        #three items
        input1 = [abc, abc, acb]
        self.assertEqual(tb.get_common_lineage(input1), ['a'])
        self.assertEqual(tb.get_common_lineage(input1, threshold=0.6), abc)
        #single item
        self.assertEqual(tb.get_common_lineage([abc]), abc)
        self.assertEqual(tb.get_common_lineage([abc], threshold=0.6), abc)
        #empty
        self.assertEqual(tb.get_common_lineage([[]]), [])
        self.assertEqual(tb.get_common_lineage([[]], threshold=0.6), [])


class FunctionTests(TestCase):
    Verbose = True
    def setUp(self):
        self.abc = ['a', 'ab', 'abc']
        self.acb =  ['a', 'ac', 'acb']
        pass
    def test_each_query_gis(self):
        obs = list(each_query_gis(StringIO(BLAST_3)))
        if self.Verbose: self.p(obs)
        self.assertEqual(len(obs), 3)
        self.assertEqual(len(obs[0][1]), 3)
        self.assertEqualItems(obs[0][1], ['58566040', '89112992', '37777029'])
        self.assertEqualItems(obs[0][0], 'G910P31FA1.T0')

    def test_each_gi_taxon_id(self):
        result = list(each_gi_taxonId(['60133320', '60133321', '60133322']))
        self.p(result)
        self.assertEqual(len(result), 3)

    #need online
    def _test_each_taxid_lineage(self):
        expect = {'10090': ['cellular organisms', 'Eukaryota',
            'Fungi/Metazoa group', 'Metazoa', 'Eumetazoa', 'Bilateria',
            'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata',
            'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii',
            'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria',
            'Euarchontoglires', 'Glires', 'Rodentia', 'Sciurognathi',
            'Muroidea', 'Muridae', 'Murinae', 'Mus'],
            '9606': ['cellular organisms', 'Eukaryota', 'Fungi/Metazoa group',
            'Metazoa', 'Eumetazoa', 'Bilateria', 'Coelomata',
            'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata',
            'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii',
            'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria',
            'Euarchontoglires', 'Primates', 'Haplorrhini', 'Simiiformes',
            'Catarrhini', 'Hominoidea', 'Hominidae', 'Homo/Pan/Gorilla group',
            'Homo']}
        result = dict(each_taxid_lineage(['10090', '9606']))
        self.assertEqual(result, expect)

    def test_each_seqid_lineage_from_european_rdp(self):
        result = each_seqid_lineage_from_european_rdp([StringIO(EMBL_2)])
        expect = [('Pla._environmental_clone_OCS162_AF001659',
            ['Plastid', 'environmental samples']), 
            ('Pla._Lygodium_japonicum_U24588',
            ['Plastid', 'Viridiplantae', 'Streptophyta', 'Embryophyta',
                'Tracheophyta', 'Filicophyta', 'Filicopsida',
                'Filicales', 'Schizaeaceae', 'Lygodium'])]
        self.assertEqual(list(result), expect)

    def test_each_query_lineages(self):
        abc, acb = self.abc, self.acb
        eid_lineages = {\
                'Fusarium_cerealis_AF141947': abc,
                'Cordyceps_prolifica_AB027324':acb,
                }
        result_dict = dict(each_query_lineages(
                [StringIO(BLAST_E_2)], eid_lineages))

        self.assertEqual(len(result_dict), 2)
        result_1 = result_dict['G912P34FA1.T0']
        self.assertEqual(result_1, #'Fusarium_cerealis_AF141947'],
                [abc, acb, acb,  abc])
        #self.assertEqual(result_1[1], #'Cordyceps_prolifica_AB027324'],
        #        ['a', 'ac', 'acb'])
        #self.assertEqual(result_1[4], #'Glomerobolus_gelineus_AF050097'],
        #        None)

    def test_each_query_toplineages(self):
        abc, acb = self.abc, self.acb
        eid_lineages = {\
                'Fusarium_cerealis_AF141947': abc,
                'Cordyceps_prolifica_AB027324':acb,
                }
        result_dict = dict(each_query_toplineages(
                StringIO(BLAST_E_2), eid_lineages))

        self.assertEqual(len(result_dict), 2)
        
        self.assertEqual(result_dict['G912P34FA1.T0'],
                [abc, acb, acb,  abc])
        self.assertEqual(result_dict['G912P34RJ21.T0'],
                [acb])


BLAST_E_2 = """\
# BLASTN 2.2.15 [Oct-15-2006]
# Query: G912P34FA1.T0
# Database: all_ssu
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
G912P34FA1.T0	Fusarium_cerealis_AF141947	96.18	838	21	11	42	877	342	1170	0.0	1320
G912P34FA1.T0	Cordyceps_prolifica_AB027324	96.18	838	21	11	42	877	521	1349	0.0	1320
G912P34FA1.T0	Cordyceps_prolifica_AB027324	96.18	838	21	11	42	877	521	1349	0.0	1320
G912P34FA1.T0	Fusarium_cerealis_AF141947	96.18	838	21	11	42	877	342	1170	0.0	1320
# BLASTN 2.2.15 [Oct-15-2006]
# Query: G912P34RJ21.T0
# Database: all_ssu
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
G912P34FA1.T0	Cordyceps_prolifica_AB027324	96.18	838	21	11	42	877	521	1349	0.0	1320
G912P34FA1.T0	Fusarium_cerealis_AF141947	96.18	838	21	11	42	877	342 1170	0.0	1120
"""
EMBL_2 = """\
ID   Pla._environmental clone OCS162  AF001659
XX
AC   AF001659; 
XX
OS   Pla. environmental clone OCS162 AF001659
OC   Plastid; environmental samples
XX
RN   [1] SSU
RA   Rappe M.S., Suzuki M.T., Vergin K.L., Giovannoni S.J.;
RT   "Phylogenetic diversity of ultraplankton plastid small-subunit rRNA genes recovered in environmental nucleic acid samples from the Pacific and Atlantic coasts of the United States";
RL   Appl. Environ. Microbiol. 64:294-303(1998).
XX
FH   Key             Location/Qualifiers
FT   source          1..6561
FT                   /strain=""
FT                   /source=""
FT                   /organism="Pla. environmental clone OCS162 AF001659"
XX
SQ  Sequence 6561 BP;
    ---------- ---------- ---------- ---------- ---------- ----------
    ---------- ---------- ---------- ---------- ------ooAG A-GUUU----
    -GAU-CAU-G
//
ID   Pla._Lygodium japonicum U24588
XX
AC   U24588; 
XX
OS   Pla. Lygodium japonicum U24588
OC   Plastid; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Filicophyta; Filicopsida; Filicales; Schizaeaceae; Lygodium
XX
RN   [1] SSU
RA   Manhart J.R.;
RT   "Chloroplast 16S rDNA sequences and phylogenetic relationships of  fern allies and ferns";
RL   Unpublished :().
XX
FH   Key             Location/Qualifiers
FH
FT   rRNA            1..6561
FT                   /product="SSU"
FT                   /accession="U24588"
FT                   /change=""
FT                   /remark=""
FT   source          1..6561
FT                   /strain=""
FT                   /source=""
FT                   /organism="Pla. Lygodium japonicum U24588"
XX
SQ  Sequence 6561 BP;
    ---------- ---------- ---------- ---------- ---------- ----------
    ---------- ---------- ---------- ---------- ------oooo o-oooo----
    -ooo-ooo-o -oooo-oooo -oooo--ooC U---GGCG-G CA-U-GC--- -C-UAA--CA
    -CA--UGC-A --AGUC-GG- -ACGG-GGAG -CACC----- ---------- ----------
//
"""

BLAST_3 = """\
# BLASTN 2.2.15 [Oct-15-2006]
# Query: G910P31FA1.T0
# Database: nt
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
G910P31FA1.T0	gi|89112992|gb|DQ401104.1|	99.26	544	2	2	38	579	511	1054	0.0	1045
G910P31FA1.T0	gi|89112992|gb|DQ401104.1|	99.59	242	0	1	593	833	1067	1308	4e-129	 470
G910P31FA1.T0	gi|37777029|dbj|AB108787.1|	99.26	544	2	2	38	579	521	1064	0.0	1045
G910P31FA1.T0	gi|37777029|dbj|AB108787.1|	99.59	242	0	1	593	833	1077	1318	4e-129	 470
G910P31FA1.T0	gi|58566040|gb|AY887946.1|	99.26	544	2	2	38	579	159	702	0.0	1045
# BLASTN 2.2.15 [Oct-15-2006]
# Query: G910P31RA1.T0
# Database: nt
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
G910P31RA1.T0	gi|45581383|dbj|AB167384.1|	88.79	232	17	9	73	297	1602	1373	4e-61	 244
G910P31RA1.T0	gi|37777029|dbj|AB108787.1|	88.79	232	17	9	73	297	1582	1353	4e-61	 244
G910P31RA1.T0	gi|59895226|gb|AY897965.1|	88.79	232	17	9	73	297	23	252	4e-61	 244
G910P31RA1.T0	gi|117653160|gb|EF060755.1|	90.86	186	10	7	73	253	186	3	1e-55	 226
G910P31RA1.T0	gi|117653104|gb|EF060699.1|	91.06	179	9	7	73	246	186	10	6e-54	 220
G910P31RA1.T0	gi|58566040|gb|AY887946.1|	88.21	212	17	8	92	297	1200	991	4e-52	 214
G910P31RA1.T0	gi|89112992|gb|DQ401104.1|	88.10	210	17	8	94	297	1550	1343	6e-51	 210
G910P31RA1.T0	gi|116174501|dbj|AB277860.1|	87.69	195	12	10	73	258	1964	1773	3e-43	 184
G910P31RA1.T0	gi|109638190|dbj|AB263744.1|	87.69	195	12	10	73	258	3474	3283	3e-43	 184
G910P31RA1.T0	gi|7643833|gb|AF150492.1|AF150492	87.69	195	12	10	73	258	882	691	3e-43	 184
# BLASTN 2.2.15 [Oct-15-2006]
# Query: G910P31FA2.T0
# Database: nt
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
G910P31FA2.T0	gi|18032809|gb|AF316593.1|	99.30	854	2	4	43	895	492	1342	0.0	1645
G910P31FA2.T0	gi|83778434|gb|DQ316826.1|	99.30	854	2	4	43	895	513	1363	0.0	1641
G910P31FA2.T0	gi|80973653|gb|DQ232535.1|	99.30	854	2	4	43	895	487	1337	0.0	1641
G910P31FA2.T0	gi|80973651|gb|DQ232533.1|	99.30	854	2	4	43	895	485	1335	0.0	1641
G910P31FA2.T0	gi|83339762|gb|DQ303190.1|	99.30	854	2	4	43	895	496	1346	0.0	1641
G910P31FA2.T0	gi|39985507|gb|AY485602.1|	99.30	854	2	4	43	895	495	1345	0.0	1641
G910P31FA2.T0	gi|30908834|gb|AY278632.1|	99.30	854	2	4	43	895	493	1343	0.0	1641
G910P31FA2.T0	gi|30908811|gb|AY278609.1|	99.30	854	2	4	43	895	495	1345	0.0	1641
G910P31FA2.T0	gi|32396627|gb|AY281086.1|	99.30	854	2	4	43	895	495	1345	0.0	1641
G910P31FA2.T0	gi|23266568|gb|AY134908.1|	99.30	854	2	4	43	895	496	1346	0.0	1641
"""
if __name__ == '__main__':
    main()

