table postarData
"RBP binding sites on the gene"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"Name of gene"
uint    score;		"Score"
char[1] strand;		"+ or - for strand"
uint    thickStart;	"Coding region start"
uint    thickEnd;	"Coding region end"
uint  	reserved;	"Green on + strand, Red on - strand"
lstring	postarID;	"POSTAR database ID"
lstring	dataSource;	"Data Source"
lstring	cellType;	"Cell type"
lstring	expSource;	"experimental source"
lstring	postarScore;	"score"
)