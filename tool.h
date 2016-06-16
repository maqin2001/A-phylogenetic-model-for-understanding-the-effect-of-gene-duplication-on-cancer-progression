extern char		spacer[10];
extern char		*printString;           
extern size_t	printStringSize;           

/* tool functions*/
FILE	*gfopen(char *filename, char *mode);
void	SetSeed (unsigned int seed);
double	LnGamma (double x);
double	rndu (void);
void	MrBayesPrint (char *format, ...);
void 	*SafeMalloc(size_t s);
double	lnchoose (int n, int x);

/* tree functions */
int			ReadaTree (FILE *fTree,Tree *tree);
void		PrintTreeToFile (FILE *file, Tree *tree);
int			PrintTree (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted);
void 		PrintNodeToFile (FILE *file, Tree *tree);
void 		findOffsprings (int *offsprings, Tree *tree, int inode);
double 		TreeHeight (Tree *tree);
int 		addNode (Tree *tree, int fromnode, int fromfather, int tonode);
int 		deleteNode (Tree *tree, int inode);
int 		swapNodes (Tree *tree, int inode, int jnode);
int 		rearrangeNodes (Tree *tree, int fromnode, int tonode);
void 		copyTree(Tree *from, Tree *to);