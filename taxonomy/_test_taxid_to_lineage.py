### test data
nodes_test = """\
1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
2	|	1	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
6	|	335928	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
7	|	6	|	species	|	AC	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
9	|	32199	|	species	|	BA	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
10	|	2	|	no rank	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
12	|	10	|	species	|	CG	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
13	|	203488	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
14	|	13	|	species	|	DT	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
16	|	32011	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
17	|	16	|	species	|	MM	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
18	|	10	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
19	|	18	|	species	|	PC	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
20	|	76892	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
""".splitlines()
names_test = """\
1	|	all	|		|	synonym	|
1	|	root	|		|	scientific name	|
2	|	Bacteria	|	Bacteria <bacteria>	|	scientific name	|
2	|	Monera	|	Monera <Bacteria>	|	in-part	|
2	|	Procaryotae	|	Procaryotae <#1>	|	in-part	|
2	|	Prokaryota	|	Prokaryota <#2>	|	in-part	|
2	|	Prokaryotae	|	Prokaryotae <#1>	|	in-part	|
2	|	eubacteria	|		|	genbank common name	|
2	|	eubacteria	|	Bacteria<blast2>	|	blast name	|
2	|	not Bacteria Haeckel 1894	|		|	synonym	|
6	|	Azorhizobium	|		|	scientific name	|
6	|	Azorhizobium Dreyfus et al. 1988	|		|	synonym	|
6	|	Azotirhizobium	|		|	equivalent name	|
7	|	Azorhizobium caulinodans	|		|	scientific name	|
7	|	Azorhizobium caulinodans Dreyfus et al. 1988	|		|	synonym	|
7	|	Azotirhizobium caulinodans	|		|	equivalent name	|
9	|	Acyrthosiphon pisum symbiont P	|		|	includes	|
9	|	Buchnera aphidicola	|		|	scientific name	|
9	|	Buchnera aphidicola Munson et al. 1991	|		|	synonym	|
10	|	"Cellvibrio" Winogradsky 1929	|		|	synonym	|
10	|	Cellvibrio	|		|	scientific name	|
10	|	Cellvibrio (ex Winogradsky 1929) Blackall et al. 1986 emend. Humphry et al. 2003	|		|	synonym	|
11	|	'Cellvibrio gilvus'	|		|	synonym	|
11	|	Cellvibrio gilvus	|		|	scientific name	|
13	|	Dictyoglomus	|		|	scientific name	|
13	|	Dictyoglomus Saiki et al. 1985	|		|	synonym	|
14	|	Dictyoglomus thermophilum	|		|	scientific name	|
14	|	Dictyoglomus thermophilum Saiki et al. 1985	|		|	synonym	|
16	|	Methyliphilus	|		|	equivalent name	|
16	|	Methylophilus	|		|	scientific name	|
16	|	Methylophilus Jenkins et al. 1987	|		|	synonym	|
16	|	Methylotrophus	|		|	misspelling	|
17	|	Methyliphilus methylitrophus	|		|	equivalent name	|
17	|	Methyliphilus methylotrophus	|		|	equivalent name	|
17	|	Methylophilus methylitrophus	|		|	equivalent name	|
17	|	Methylophilus methylotrophus	|		|	scientific name	|
17	|	Methylophilus methylotrophus Jenkins et al. 1987	|		|	synonym	|
17	|	Methylotrophus methylophilus	|		|	synonym	|
18	|	Pelobacter	|		|	scientific name	|
18	|	Pelobacter Schink and Pfennig 1983	|		|	synonym	|
19	|	Pelobacter carbinolicus	|		|	scientific name	|
19	|	Pelobacter carbinolicus Schink 1984	|		|	synonym	|
20	|	Phenylobacterium	|		|	scientific name	|
20	|	Phenylobacterium Lingens et al. 1985 emend. Kanso and Patel 2004	|		|	synonym	|
20	|	Phenylobacterium Lingens et al. 1985 emend. Tiago et al. 2005	|		|	synonym	|
""".splitlines()

