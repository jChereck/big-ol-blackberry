#include "mat.h"
#include <float.h>

// Uses a Matrix to create a balanced binary tree from rows of data
// by repeatedly sorting on different columns.   The first column
// is the label of the row which is a known classification or property of the
// data in the rest of the row.   This fits well with the labeled row version
// supported by the matrix library.
//
// tree: tree to be built in Matrix form
// rowstart, rowend: the range of rows to sort
// c: the column to sort on
void build(Matrix &tree, int rowstart, int rowend, int c)
{
    int split;

    if (c>=tree.numCols()) c=1;

    tree.sortRowsByCol(c, rowstart, rowend);

    //printf("start: %d end: %d c: %d\n", rowstart, rowend, c);
    //tree.print();

    split = (rowstart + rowend)/2;
    if ((split-1)-rowstart>0) build(tree, rowstart, split-1, c+1);
    if (rowend-(split+1)>0) build(tree, split+1, rowend, c+1);
}



// IMPORTANT: must be full real distance measure (not the square of
// the distance) since short cut relies on max abs(difference) as a
// distance in one dimension. ASSUMES matrix is a labeled matrix.
double myDist(Matrix &x, int r, Matrix &item) {
	//x.subMatrix(r,1,1,0).print();
	//item.print();
    return sqrt(x.subMatrix(r, 1, 1, 0).dist2(item));
}

int compElm(Matrix &tree, int r, Matrix &item, int c){
	
	double it = item.get(0,c-1);
	double tr = tree.get(r,c);

	if       ( it == tr ){
		//printf("%f = %f\n", it, tr);
		return 0;
	}else if ( it < tr ){
		//printf("%f < %f\n", it, tr);
		return -1;
	}else{
		//printf("%f > %f\n", it, tr);
		return 1;
	}
}

// tree is a kd-tree and item is a row matrix to look up.
void nearestAux(Matrix &tree,   // the kdtree matrix
                Matrix &item,   // the row
                int rowstart,   // beginning row of region to search (subtree)
                int rowend,     // ending row of region to search
                int c,          // column or feature to compare
                double &best,   // the distance to nearest point found so far
                int &bestrow) {  // the row of the nearest point found so far

	//if we are at the base of the tree, brute force
	if( c >= tree.numCols()-1 ){
		//printf("Done with tree traverse\n");
		for( int r = rowstart; r <= rowend; r++ ){
			double dist = myDist(tree, r, item);
			//printf("dist is %f\n", dist);
			if( dist < best ){
				//printf("new shortest!\n");
				best = dist;
				bestrow = r;
			}
		}

		//one you find the best row, return
		return;
	}

	//where is split for tree
	int split = (rowstart + rowend)/2;
	
	//test current split point
	double dist = myDist(tree, split, item);
	//printf("dist is %f\n", dist);
	if( dist < best ){
		//printf("new shortest!\n");
		best = dist;
		bestrow = split;
	}

	//recurse down proper path
	int comp = compElm(tree, split, item, c);
	if( comp >= 0 ){
		//printf("recursing down from split: %d\n", split);
		nearestAux(tree, item, split+1, rowend, c+1, best, bestrow);
	}

	if (comp <= 0){
		//printf("recursing up from split: %d\n", split);
		nearestAux(tree, item, rowstart, split-1, c+1, best, bestrow);
	}
	

    return;
}


int nearest(Matrix &tree, Matrix &item) {
    double best;
    int bestrow;

    best = DBL_MAX;  // minimum value found so far (maximum double value from float.h)
    bestrow = -1;    // row of best
    nearestAux(tree, item, 0, tree.numRows()-1, 1, best, bestrow);

    return bestrow;
}



int main()
{
    Matrix items("item");      // the data item to be searched for
    Matrix tree("kdtree");     // a kdtree as a matrix
    SymbolNumMap labels;       // labels for a row labeled matrix (see the matrix library)
    int closest;               // the row that is closest to item using the measure myDist

    // read in labeled training data
    labels = tree.readLabeledRow();

    // build tree
    build(tree, 0, tree.numRows()-1, 1);

    // display tree for comparison in tests
    printf("KDTree version of matrix");
    tree.printLabeledRow(labels);
    printf("\n");

    // read in the items we wish to look up
    items.read();

    MatrixRowIter a(&items);
    for (Matrix *item = a.rowBegin(); a.rowNotEnd(); a.rowNext()) {

        // display the item
        printf("\nSOLVE: ");
        item->writeLine(0);
        printf("\n");

        // use kdtree to find row number of closest in data
        closest = nearest(tree, *item);

        printf("ANS: ");
        printf("%d %s ", closest, labels.getStr(tree.get(closest, 0)).c_str());
        tree.subMatrix(closest, 1, 1, 0).writeLine(0);
        printf("\n");
    }

    return 0;
}

