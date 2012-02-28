/****************************************************************

This Slice_iter is based on that defined in Stroustrup 

                   A. Bartholomew 3rd January, 2005

*****************************************************************/


namespace Slice_iterators {
  using std::valarray;
  using std::slice;

  template<class T> class Slice_iter{
    valarray<T>* v;
    slice s;
    size_t curr;   //index of current element
	  
    T& ref(size_t i) const {return (*v)[s.start()+i*s.stride()];}
	
  public:
    Slice_iter(valarray<T>* vv, slice ss): v(vv), s(ss), curr(0){}
	
	size_t size() const {return s.size();}
		
    Slice_iter end() const{
      Slice_iter t = *this;
      t.curr=s.size();
      return t;
    }
	
    Slice_iter& operator++(){curr++; return *this;}
    Slice_iter  operator++(int){Slice_iter t=*this; curr++; return t;}
    Slice_iter& operator+=(int n){curr+=n; return *this;}

    Slice_iter& operator--(){curr--; return *this;}
    Slice_iter  operator--(int){Slice_iter t=*this; curr--; return t;}
    Slice_iter& operator-=(int n){curr-=n; return *this;}

    
    T& operator[] (size_t i) const {return ref(i); }    //C style subscript  
//    T& operator() (size_t i) {return ref(i); }    //Fortran style subscript  
    T& operator*()           {return ref(curr); } //current element

    friend
      bool operator==(const Slice_iter<T>& p, const Slice_iter<T>& q){
      //      bool operator==<>(const Slice_iter<T>& p, const Slice_iter<T>& q){
      return 
        p.curr==q.curr && 
        p.s.stride()==q.s.stride() && 
        p.s.start() == q.s.start();
    }
    
    friend
      bool operator!=(const Slice_iter<T>& p, const Slice_iter<T>& q){
      //      bool operator!=<>(const Slice_iter<T>& p, const Slice_iter<T>& q){
      return  (!p==q);
    }
    
    friend
      bool operator< (const Slice_iter<T>& p, const Slice_iter<T>& q){
      //      bool operator< <>(const Slice_iter<T>& p, const Slice_iter<T>& q){
      return 
        p.curr<q.curr && 
        p.s.stride() == q.s.stride() && 
        p.s.start()  == q.s.start();
    }
  };

  template<class T> class Cslice_iter{
    valarray<T>* v;
    slice s;
    size_t curr;   //index of current element
    T& ref(size_t i) const {return (*v)[s.start()+i*s.stride()];}
  public:
    Cslice_iter(valarray<T>* vv, slice ss): v(vv), s(ss), curr(0){}
    Cslice_iter end() const{
      Cslice_iter t = *this;
      t.curr=s.size();
      return t;
    }
    const Cslice_iter& operator++(){curr++; return *this;}
    Cslice_iter  operator++(int){Cslice_iter t=*this; curr++; return t;}
    
    const T& operator[] (size_t i) const  {return ref(i); }//C style subscript  
//    const T& operator() (size_t i) const {return ref(i); }    //Fortran style subsc  
    const T& operator*()           const {return ref(curr); } //current element

    friend
      bool operator == (const Cslice_iter<T>& p, const Cslice_iter<T>& q){
      //      bool operator ==<> (const Cslice_iter<T>& p, const Cslice_iter<T>& q){
      return 
        p.curr==q.curr && 
        p.s.stride()==q.s.stride() && 
        p.s.start() == q.s.start();
    }
    
    friend 
      //      bool operator !=<> (const Cslice_iter& p, const Cslice_iter& q){
      bool operator != (const Cslice_iter& p, const Cslice_iter& q){
      return  (!p==q);
    }
 
    friend
      bool operator < (const Cslice_iter& p, const Cslice_iter& q){
      //      bool operator < <>(const Cslice_iter& p, const Cslice_iter& q){
      return 
        p.curr<q.curr && 
        p.s.stride() == q.s.stride() && 
        p.s.start()  == q.s.start();
    }
    
  };

  template<class T> void print (Slice_iter<T> r, ostream& s, int n, string prefix)
  {
	  s << prefix;
	  for (size_t i=0; i< r.size(); i++)
	  {
		  if (n !=0)
			  s << setw(n) << r[i];
		  else
			  s << r[i] << ' ';
	  }
	  s << endl;
  }

}//namespace

