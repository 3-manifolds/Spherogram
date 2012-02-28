/* the Scalar type St is the type of the scalar necessary to initialize a T with
   zero, so we can write e.g. T det = T(St(0));  This is necessary to overcome the
   restriction of one implicit conversion between user defined types.
*/
template<class T, class St = int> class matrix
{

    valarray<T>* v;
	size_t rows;
	size_t cols;

public:
	/* These are friends to allow implicit conversion of their arguments,
       helper functions do not support implicit conversion */
	friend matrix<T,St> operator + <> (const matrix<T,St>&, const matrix<T,St>&);
	friend matrix<T,St> operator - <> (const matrix<T,St>&, const matrix<T,St>&);
	friend matrix<T,St> operator * <> (const matrix<T,St>&, const matrix<T,St>&);
	friend matrix<T,St> operator * <> (const T&, const matrix<T,St>&);
	friend bool operator != <> (const matrix<T,St>&, const matrix<T,St>&);
	friend ostream& operator << <> (ostream& outstream, const matrix<T,St>&);
	
	friend matrix<T,St> operator - <> (const int&, const matrix<T,St>&);
	friend bool operator == <> (const matrix<T,St>& , const int );
	friend bool operator != <> (const matrix<T,St>&, const int );
	friend bool operator == <> (const int, const matrix<T,St>&);
	friend bool operator != <> (const int, const matrix<T,St>&);

	friend T trace <> (const matrix<T,St> M);
	friend T determinant <> (const matrix<T,St>& M, string title, int n, int* rperm, int* cperm, int recursion_level);
	
	typedef St scalar_type;
    matrix<T,St> (size_t r, size_t c);
	matrix<T,St> ( const matrix<T,St>&); // copy constructor
    ~matrix<T,St>() {delete v;}
	size_t size() const {return rows*cols;}
    size_t numrows() const {return rows;}
    size_t numcols() const {return cols;}
	valarray<T>& array() {return *v;}

	Slice_iter<T> row (size_t i);
	Cslice_iter<T> row (size_t i) const;
	Slice_iter<T> column (size_t i);
	Cslice_iter<T> column (size_t i) const;

	Slice_iter<T> operator[] (size_t i) {return row(i);}
	Cslice_iter<T> operator[] (size_t i) const {return row(i);}

    matrix<T,St>& operator = (const matrix<T,St>&); // copy assignment
	matrix<T,St> operator += (const matrix<T,St>&);
	matrix<T,St> operator -= (const matrix<T,St>&);
	matrix<T,St> operator *= (const matrix<T,St>&);
	matrix<T,St> operator *= (const T&);
	matrix<T,St> inverse(bool adjunct_only) const;
	bool operator == (const matrix<T,St>&) const;
	matrix<T,St> operator -= (const int);
	void dump(ostream& os) const;
};

template <class T, class St> 
void set_matrix_N_element(matrix<T,St>& d_matrix, int d_N_row, int d_N_col, 
                          const matrix<T,St>& s_matrix, int s_N_row, int s_N_col, int N);

template <class T, class St> 
void decrement_matrix_N_element(matrix<T,St>& d_matrix, int d_N_row, int d_N_col, int N);
                          
template <class T, class St> 
void print(matrix<T,St>& m, ostream& s, int n, string prefix);
