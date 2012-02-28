/**********************************************************************
  matrix.h contains a matrix template definition based on the Matrix
  class described by Strustrup.

  NOTE: This representation of matrices stores them row by row in a valarray
        rather than column by column in a valarray.  This is because I think
		a row by row approach is more natural 

  this matrix.h was derived from the array-of-pointers template version
            A. Bartholomew, 3rd January, 2005
			
***********************************************************************/
#include <valarray>
#include <Slice_iter.h>

using namespace Slice_iterators;

struct matrix_control 
{
	static unsigned int DEBUG; //bitmap
	enum parameters
	{
		general=1, // unused
		determinant=2, 
		multiply=4, 
		inverse=8, 
		// 'all' also supported
	};

	static bool WAIT_INFO;
	static bool COMFORT_DOTS;
	static bool SINGLE_LINE_OUTPUT;
	static int DET_DEBUG_LIMIT;
	static int wait_threshold;
	static int wait_count;
	static int reset_count; // number of times wait_count has reached wait_threshold
};

struct matrix_error {
	matrix_error (string message) {cout << "\nMatrix error!" << message << endl;}
};


/* the Scalar type St is the type of the scalar necessary to initialize a T with
   zero, so we can write e.g. T det = T(St(0));  This is necessary to overcome the
   restriction of one implicit conversion between user defined types.
*/
template<class T, class St> class matrix;

/* First declare the friend functions */

template<class T, class St> inline matrix<T,St> operator + (const matrix<T,St>& a, const matrix<T,St>& b)
{
	matrix<T,St> result = a;
	return result += b;
}

template<class T, class St> inline matrix<T,St> operator - (const matrix<T,St>& a, const matrix<T,St>& b)
{
	matrix<T,St> result = a;
	return result -= b;
}

template<class T, class St> inline matrix<T,St> operator * (const matrix<T,St>& a, const matrix<T,St>& b)
{
	matrix<T,St> result = a;
	return result *= b;
}

template<class T, class St> inline matrix<T,St> operator * (const T& t, const matrix<T,St>& b)
{
	matrix<T,St> result = b;
	return result *= t;
}

template<class T, class St> inline bool operator != (const matrix<T,St>& a, const matrix<T,St>& b)
{
	return !(a==b);
}

template<class T, class St> ostream& operator << (ostream& outstream, const matrix<T,St>& m)
{
    if (!matrix_control::SINGLE_LINE_OUTPUT)
		outstream << "\n";
		
    for (size_t i = 0; i< m.numrows(); i++ )
    {
        for (size_t j = 0; j < m.numcols(); j++)
			outstream << m[i][j] << ' ';

	    if (!matrix_control::SINGLE_LINE_OUTPUT)
		    outstream << endl;
    }
    return outstream;
}

/* subtract this from n x I */
template<class T, class St> matrix<T,St> operator - (const int& num, const matrix<T,St>& M)
{
	matrix<T,St> result(M.numrows(),M.numcols());
	St minus_one = -1;
	for (size_t i = 0; i< M.numrows(); i++)
	{
		for (size_t j=0; j< M.numcols(); j++)
		{
			result[i][j] = M[i][j]*T(minus_one);
		}
		
		St Mnum = num;
		result[i][i] += T(Mnum);
	}
	return result;
}

template<class T, class St> bool operator == (const matrix<T,St>& M, const int num )
{
	for (int i = 0; i< M.numrows(); i++)
	{
		for (int j=0; j< i; j++)
		{
			if (M[i][j] != T(St(0)))
				return false;
		}
		if (M[i][i] != T(St(num)))
			return false;
		for (int j=i+1; j< M.numcols(); j++)
		{
			if (M[i][j] != T(St(0)))
				return false;
		}
	}
	return true;
}

template<class T, class St> inline bool operator != (const matrix<T,St>& M, const int num )
{
	return !(M == num);
}

template<class T, class St> inline bool operator == (const int num, const matrix<T,St>& M)
{
	return M == num;
}

template<class T, class St> inline bool operator != (const int num, const matrix<T,St>& M)
{
	return !(M == num);
}

template<class T, class St> inline T trace (const matrix<T,St> M)
{
	T result;
	
	for (int i = 0; i < M.numrows(); i++)
		result += M[i][i];
	
	return result;
}

/* determinant uses a row and column permutation rperm and cperm, and the number of rows/columns n, 
   to evaluate the determinant of any sub-matrix of the square matrix M.  The recursion_level is used
   only for debugging purposes.

   The approach is to look along the top row and down the left column of the submatrix determined by rperm 
   and cperm, and to evaluate the determinant based on the smallest number of non-zero entries.  The function 
   then recurses to evaluate the sub-determinants; only when we get to n = 2 do we have a 2x2 submatrix 
   whereupon we calculate the determinant directly.
*/

template <class T, class St> 
T determinant (const matrix<T,St>& M, string title="untitled", int n=0, int* rperm=0, int* cperm=0, int recursion_level=0)
{

	// last change 3/9/06
	
    T zero = T(St(0));
    T det = zero;
    T temp = zero;
	bool clean_up_perms = false;
	

	/* We are starting a new determinant calculation if recursion_level == 0 */
	if (!recursion_level)
	{
		matrix_control::reset_count = 0;
if (matrix_control::DEBUG & matrix_control::determinant)
	debug << "matrix::determinant: underlying matrix M = \n" << M << endl;
	}
	
	if (rperm == 0)
	{
	
if (matrix_control::DEBUG & matrix_control::determinant)
	debug << "matrix::determinant: default parameters provided, creating permutations" << endl;
	
		n = M.numcols();
		rperm = new int[M.numrows()];
		for (size_t i=0; i< M.numrows(); i++)
			rperm[i] = i;
		cperm = new int[n];
		for (int i=0; i< n; i++)
			cperm[i] = i;
		clean_up_perms = true;
	}


if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "n = " << n << ", recursion_level = " << recursion_level << endl;

    debug << "matrix::determinant: ";
	for (int k=0;k<recursion_level;k++)
	    debug << "    ";
	debug << "rperm: ";
	for ( int k = 0 ; k< n; k++)
    	debug << rperm[k] << " ";
	debug << endl;

    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "cperm: ";
    for (int k = 0 ; k< n; k++)
		debug << cperm[k] << " ";
	debug << endl;
}

	if (n==1)
	{

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
	debug << "matrix::determinant: returning single element " << M[rperm[0]][cperm[0]] << endl;

		if (clean_up_perms)
		{
			delete[] rperm;
			delete[] cperm;
		}
		return M[rperm[0]][cperm[0]];
	}
    else if (n == 2)
    {
if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
	for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "positive term: [" << rperm[0] << "][" << cperm[0] << "]*["
          << rperm[1] << "][" << cperm[1] << "] = " 
	      << M[rperm[0]][cperm[0]] * M[rperm[1]][cperm[1]];
	debug << endl;
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "negative term: [" << rperm[0] << "][" << cperm[1] << "]*["
          << rperm[1] << "][" << cperm[0] << "] = "
	      << M[rperm[0]][cperm[1]] * M[rperm[1]][cperm[0]];
	debug << endl;	
}
		det =   M[rperm[0]][cperm[0]] * M[rperm[1]][cperm[1]]
		      - M[rperm[0]][cperm[1]] * M[rperm[1]][cperm[0]];

    }
    else
    {
		int sub_r_perm[n-1];
		int sub_c_perm[n-1];
		bool evaluate_along_row = true;
		
		/* Look to see whether there are more zeros along the top row or left column
		   we evaluate the determinant based on the largest number of zero values
		*/
		int num_zeros = 0;
		for (int i=0; i < n; i++)
		{
			if (M[rperm[0]][cperm[i]] == zero)
				num_zeros++;
				
			if (M[rperm[i]][cperm[0]] == zero)
				num_zeros--;
		}
				
		if (num_zeros < 0)
			evaluate_along_row = false;

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "evaluating " << n << "-sub-determinant " << (evaluate_along_row? "along row" : "down column") 
	      << ", row-column num_zeros = " << num_zeros << endl;
}

		if (evaluate_along_row) // row permutation for all recursive calls at this level is the same
		{
			for (int i=1; i<n;i++)
	    		sub_r_perm[i-1] = rperm[i];

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "sub_r_perm (fixed for iterative calls at this level): ";
	for ( int k = 0 ; k< n-1; k++)
	   	debug << sub_r_perm[k] << " ";
	debug << endl;
}			
		}
		else // column permutation for all recursive calls at this level is the same
		{
			for (int i=1; i<n;i++)
	    		sub_c_perm[i-1] = cperm[i];
			
if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "sub_c_perm (fixed for iterative calls at this level): ";
	for ( int k = 0 ; k< n-1; k++)
	   	debug << sub_c_perm[k] << " ";
	debug << endl;
}			
		}
				
		for (int i=0; i<n; i++)
		{
		
		    if ((evaluate_along_row && M[rperm[0]][cperm[i]] != T(St(0))) || (!evaluate_along_row && M[rperm[i]][cperm[0]] != T(St(0))))
	    	{
			
				if (evaluate_along_row)
				{
					for (int j=0; j<i;j++)
			    		sub_c_perm[j] = cperm[j];
					for (int j=i+1;j<n;j++)
			    		sub_c_perm[j-1] = cperm[j];

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
	for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n-1 << "-multiplier =" << (i%2? " -1 * " :" ") << M[rperm[0]][cperm[i]] << endl;
}

					temp = M[rperm[0]][cperm[i]] * determinant(M,title, n-1, sub_r_perm, sub_c_perm, recursion_level+1);

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "additive term = " << temp << endl;
}

				}
				else
				{
					for (int j=0; j<i;j++)
			    		sub_r_perm[j] = rperm[j];
					for (int j=i+1;j<n;j++)
			    		sub_r_perm[j-1] = rperm[j];

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
	for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n-1 << "-multiplier =" << (i%2? " -1 * " :" ") << M[rperm[i]][cperm[0]] << endl;
}

					temp = M[rperm[i]][cperm[0]] * determinant(M,title, n-1, sub_r_perm, sub_c_perm, recursion_level+1);

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "additive term = " << temp << endl;
}
				}
				

				if (i%2)
					det -= temp;
				else
					det += temp;

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n << "-sub-determinant: " << det << endl;
}

	    	}
		}
    }

	if (matrix_control::WAIT_INFO)
	{
	    if (n == matrix_control::wait_threshold)
	    {
			if (matrix_control::COMFORT_DOTS)
			{
				cout << ".";
   				cout.flush();
			}
			
			if (matrix_control::wait_count > 1000)
			{
				matrix_control::reset_count++;
	    		cout << "\nworking on determinants for " << title << " (" << matrix_control::reset_count << "), please wait\n";
	    		matrix_control::wait_count = 0;
			}
			else
	    		matrix_control::wait_count++;
    	}
	}

if (matrix_control::DEBUG & matrix_control::determinant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::determinant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n << "-determinant: " << det << endl;
}

	if (clean_up_perms)
	{
		delete[] rperm;
		delete[] cperm;
	}
    return det;
}


#include<matrix.int.h>


template<class T, class St> matrix<T,St>::matrix (size_t r, size_t c)
{
    v = new valarray<T>(T(St(0)),r*c);
    rows = r;
    cols = c;
}

/* the copy constructor */
template<class T, class St> matrix<T,St>::matrix (const matrix<T,St>& M)
{
	cols = M.cols;
	rows = M.rows;
	
    v = new valarray<T>(T(St(0)),rows*cols);
	*v = *M.v;
}

/* row and column definitions - stored row by row in the representation */
template <class T, class St> inline Slice_iter<T> matrix<T,St>::row(size_t i)
{
	return Slice_iter<T>(v,slice(i*cols, cols, 1));
}

template <class T, class St> inline Cslice_iter<T> matrix<T,St>::row(size_t i) const
{
	return Cslice_iter<T>(v,slice(i*cols, cols, 1));
}

template <class T, class St> inline Slice_iter<T> matrix<T,St>::column(size_t i)
{
	return Slice_iter<T>(v,slice(i,rows,cols));
}

template <class T, class St> inline Cslice_iter<T> matrix<T,St>::column(size_t i) const
{
	return Cslice_iter<T>(v,slice(i,rows,cols));
}

/* copy assignment */
template<class T, class St> inline matrix<T,St>& matrix<T,St>::operator = (const matrix<T,St>& M)
{
	if (this != &M)
	{	
		cols = M.cols;
		rows = M.rows;
		
	    delete v;
		v = new valarray<T>(T(St(0)),rows*cols);
		*v = *M.v;
		
	}
	return *this;
}


/* Now the operators */
template<class T, class St> matrix<T,St> matrix<T,St>::operator += (const matrix<T,St>& M)
{
	matrix<T,St>& loc = *this;
	
	if (loc.rows != M.rows || loc.cols != M.cols)
		throw(matrix_error ("incompatible sizes in += operator"));
	
	for (size_t i=0; i< loc.numrows(); i++)
	{
		for (size_t j=0; j< loc.numcols(); j++)
		{
			loc[i][j] += M[i][j];
		}
	}
	
	return loc;
}

template<class T, class St> matrix<T,St> matrix<T,St>::operator -= (const matrix<T,St>& M)
{
	matrix<T,St>& loc = *this;
	
	if (loc.rows != M.rows || loc.cols != M.cols)
		throw(matrix_error ("incompatible sizes in -= operator"));
	
	for (size_t i=0; i< loc.numrows(); i++)
	{
		for (size_t j=0; j< loc.numcols(); j++)
		{
			loc[i][j] -= M[i][j];
		}
	}
	
	return loc;
}

template<class T, class St> matrix<T,St> matrix<T,St>::operator *= (const matrix<T,St>& M)
{
	matrix<T,St>& loc = *this;

	if (loc.rows != M.cols )
		throw(matrix_error ("incompatible sizes in -= operator"));

	matrix<T,St> result(loc.rows, M.cols);

if (matrix_control::DEBUG & matrix_control::multiply)
{
	debug << "\nmatrix::operator *= : \n\t(*this) = " << *this <<  "\n\t      M = " << M << endl;
	debug << "\nmatrix::operator *= : multiplication produces " << loc.numrows() << " by " << M.numcols() << " result";
}


	for (size_t i=0; i< loc.rows; i++)
	{
		for (size_t j=0; j< M.cols; j++)
		{
if (matrix_control::DEBUG & matrix_control::multiply)
	debug << "\nmatrix::operator *= : element " << i << "," << j <<" contributions\n";	
			for (size_t k=0; k< loc.cols; k++)
			{
				result[i][j] += loc[i][k]*M[k][j];
if (matrix_control::DEBUG & matrix_control::multiply)
{
	T temp = loc[i][k]*M[k][j];
	debug << "\t" << temp;
}
			}
		}
	}
	
	loc = result;
	return loc;
}


template<class T, class St> matrix<T,St> matrix<T,St>::operator *= (const T& t)
{
	matrix<T,St>& loc = *this;
	
	for (unsigned int i=0; i< loc.numrows(); i++)
	for (unsigned int j=0; j< loc.numcols(); j++)
		loc[i][j] *= t;

	return loc;
}


/* subtract n x I from this */
template<class T, class St> matrix<T,St> matrix<T,St>::operator -= (const int num)
{
	matrix<T,St>& loc = *this;
	matrix<T,St> result = loc;
	for (int i=0; i< loc.numrows(); i++)
	{
		result[i][i] -= T(1);
	}
	return result;
}

/* Currently this inverse function assumes that the determinant is non-zero */
template<class T, class St> matrix<T,St> matrix<T,St>::inverse(bool adjunct_only=false) const
{
	const matrix<T,St>& M = *this;
	int n = M.numrows();
	
	matrix<T,St> invM(n,n);
	
	
	T det = determinant (M);

if (matrix_control::DEBUG & matrix_control::inverse)
	debug << "\nmatrix::inverse: det(M) = " << det << endl;

	int rperm[n];
	int cperm[n];
	for (int i=0; i<n; i++)
		rperm[i]=cperm[i]=i;

	if (det == T(St(0)))
		throw(matrix_error ("attempt to take inverse of singular matrix"));
	
	/* calculate the inverse of M using the adjoint method */
	for (int i = 0; i < n; i++)
	{
		/* take out row i */
		for (int k=0; k < i; k++)
			rperm[k] = k;
		for (int k=i+1; k<n; k++)
			rperm[k-1] = k;
	
		for (int j=0; j<n; j++)
		{
			/* take out column j */
			for (int k=0; k < j; k++)
				cperm[k] = k;
			for (int k=j+1; k<n; k++)
				cperm[k-1] = k;
			
			/* element (j,i) of the adjoint is the signed n-1 x n-1 determinant of the 
			   matrix determined by rperm and cperm */
			invM[j][i] = determinant(M,"",n-1,rperm,cperm);

if (matrix_control::DEBUG & matrix_control::inverse)
	debug << "\nmatrix::inverse: initial invM[" << j << "][" << i << "] = " << invM[j][i] << endl;

			/* the sign is given by {-1}^{i+j} */
			St minus_one = St(-1);
			if ((i+j) % 2)
				invM[j][i] *= T(minus_one);

if (matrix_control::DEBUG & matrix_control::inverse)
	debug << "\nmatrix::inverse: after sign invM[" << j << "][" << i << "] = " << invM[j][i] << endl;
	
		}
	}
	
	/* divide down by the determinant to get the inverse unless we're
	   calculating the adjunct only
	*/
	if (!adjunct_only)
	{
		for (int i=0;i<n;i++)
		for (int j=0;j<n;j++)
			invM[i][j] /= det;
	}
	
	return invM;
}

template<class T, class St> bool matrix<T,St>::operator == (const matrix<T,St>& M) const
{
	const matrix<T,St>& loc = *this;
	
	if (loc.rows != M.rows || loc.cols != M.cols)
		return false;
	
	for (size_t i=0; i< loc.rows; i++)
	for (size_t j=0; j< loc.cols; j++)
	{
		if (loc[i][j] != M[i][j])
			return false;
	}
	
	return true;
}

template<class T, class St> void matrix<T,St>::dump(ostream& os) const
{
	os << "\nvalarray at " << v << " rows = " << rows << " cols = " << cols << endl;
	for (size_t r = 0; r < rows ; r++)
	{
		os << "\nrow " << r;
		for (size_t c = 0; c < cols; c++)
		{
			os << "\n\t[" << c << "] at " << &((*this)[r][c]) << " value: " << (*this)[r][c];
		}
	}
	os << endl;
}

/* The next two functions are designed for r*N x c*N matrices that may be regarded as r x c matrices having 
   N x N matrix elements.  The functions set an N x N element in a destination matrix from an N x N element
   in a source matrix, and decrement an N x N element by 1, that is subtract the N x N identity matrix
*/
template <class T, class St> 
void set_matrix_N_element(matrix<T,St>& d_matrix, int d_N_row, int d_N_col, 
                          const matrix<T,St>& s_matrix, int s_N_row, int s_N_col, int N)
{
	int d_r_base = N*d_N_row;
	int d_c_base = N*d_N_col;
	int s_r_base = N*s_N_row;
	int s_c_base = N*s_N_col;
	
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
		d_matrix[d_r_base+i][d_c_base+j] = s_matrix[s_r_base+i][s_c_base+j];
}


template <class T, class St> 
void decrement_matrix_N_element(matrix<T,St>& d_matrix, int d_N_row, int d_N_col, int N)
{
	int d_r_base = N*d_N_row;
	int d_c_base = N*d_N_col;

	for (int i = 0; i < N; i++)
		d_matrix[d_r_base+i][d_c_base+i] += T("-1");
}

template <class T, class St> 
void print(matrix<T,St>& m, ostream& s, int n, string prefix)
{
	for (unsigned int i=0; i< m.numrows(); i++)
		print (m[i],s,n,prefix);
}
