/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __math_complex_number_h__
#define __math_complex_number_h__

#include "mrtrix.h"

namespace MR {
  namespace Math {

    template <class T> class ComplexNumber {
      protected:
        T data[2];
      public:
        ComplexNumber ();
        ComplexNumber (T real, T imag);
        ComplexNumber (const ComplexNumber<T>& c);

        void                set (T real, T imag);
        T&                  re ();
        T&                  im ();

        const T&            re () const;
        const T&            im () const;
        T                   mod () const;
        T                   mod2 () const;
        T                   phase () const;

        ComplexNumber<T>    operator* (T val) const;
        ComplexNumber<T>    operator/ (T val) const;

        void                operator*= (T val);
        void                operator/= (T val);

        ComplexNumber<T>    operator+ (const ComplexNumber<T>& c) const;
        ComplexNumber<T>    operator- (const ComplexNumber<T>& c) const;
        ComplexNumber<T>    operator* (const ComplexNumber<T>& c) const;

        void                operator+= (const ComplexNumber<T>& c);
        void                operator-= (const ComplexNumber<T>& c);
        void                operator*= (const ComplexNumber<T>& c);

        ComplexNumber<T>    conj () const;
        ComplexNumber<T>    norm () const;


        ComplexNumber<T>&   operator= (const ComplexNumber<T>& c);
        T*                  pointer();
    };


    template <class T> std::ostream& operator<< (std::ostream& stream, const ComplexNumber<T>& a);









    /** \brief create instance initialised to 0 + 0\e i. 
     *
     * Default constructor: creates an instance of a complex number initialised to 0 + 0\e i. */
    template <class T> inline ComplexNumber<T>::ComplexNumber ()
    {
      set(0.0, 0.0); 
    }


    /** \brief create instance initialised to \p real + \p imag \e i. 
     *
     * Constructor: creates an instance of a complex number initialised to \p real + \p imag \e i. */
    template <class T> inline ComplexNumber<T>::ComplexNumber (T real, T imag)
    { 
      set(real, imag); 
    }


    /** \brief create instance initialised from complex number \p c. 
     *
     * Copy constructor: creates an instance of a complex number initialised from 
     * the supplied complex number \p c. */
    template <class T> inline ComplexNumber<T>::ComplexNumber (const ComplexNumber<T>& c)
    { 
      set(c.data[0], c.data[1]); 
    }


    /** \brief set complex number to \p real + \p imag \e i. 
     *
     * This function sets the complex number to \p real + \p imag \e i. */
    template <class T> inline void ComplexNumber<T>::set (T real, T imag) 
    {
      data[0] = real;
      data[1] = imag; 
    }


    /** \brief returns the real part of the complex number.
     *
     * \return a reference to the real part of the complex number.*/
    template <class T> inline T& ComplexNumber<T>::re ()
    {
      return (data[0]); 
    }


    /** \brief returns the imaginary part of the complex number.
     *
     * \return a reference to the imaginary part of the complex number. */
    template <class T> inline T& ComplexNumber<T>::im ()
    {
      return (data[1]); }


    /** \brief returns the real part of the complex number. 
     *
     * \return a const reference to the real part of the complex number. */
    template <class T> inline const T& ComplexNumber<T>::re () const
    {
      return (data[0]); 
    }


    /** \brief returns the imaginary part of the complex number. 
     * 
     * \return a const reference to the imaginary part of the complex number. */
    template <class T> inline const T& ComplexNumber<T>::im () const
    {
      return (data[1]); 
    }


    /** \brief assignment operator. 
     *
     * assign the value of the right-hand argument to the left-hand argument. */
    template <class T> inline ComplexNumber<T>& ComplexNumber<T>::operator= (const ComplexNumber<T>& c)
    {
      set (c.re(), c.im());
      return (*this); 
    }


    /** \brief returns the modulus of the complex number. 
     *
     * \return the modulus of the complex number. */
    template <class T> inline T ComplexNumber<T>::mod () const
    {
      return (sqrt(mod2())); 
    }


    /** \brief returns the squared modulus of the complex number. 
     *
     * \return the squared modulus of the complex number. */
    template <class T> inline T ComplexNumber<T>::mod2 () const
    {
      return (re()*re() + im()*im()); 
    }


    /** \brief returns the phase of the complex number. 
     *
     * \return the phase angle of the complex number, in radians. */
    template <class T> inline T ComplexNumber<T>::phase () const 
    { 
      return (atan2 (im(), re())); 
    }


    /** \brief returns a copy of the complex number multiplied by \p val. 
     *
     * \return a copy of the complex number multiplied by the real number \p val. */
    template <class T> inline ComplexNumber<T> ComplexNumber<T>::operator* (T val) const 
    {
      ComplexNumber<T> ret (val*re(), val*im()); 
      return (ret); 
    }


    /** \brief returns a copy of the complex number divided by \p val. 
     *
     * \return a copy of the complex number divided by the real number \p val. */
    template <class T> inline ComplexNumber<T> ComplexNumber<T>::operator/ (T val) const
    {
      ComplexNumber<T> ret (re()/val, im()/val); 
      return (ret);
    }


    /** \brief multiplies the complex number by \p val. 
     *
     * The complex number is multiplied by the real number \p val. */
    template <class T> inline void ComplexNumber<T>::operator*= (T val)
    {
      re() *= val;
      im() *= val; 
    }


    /** \brief divides the complex number by \p val. 
     *
     * The complex number is divided by the real number \p val. */
    template <class T> inline void ComplexNumber<T>::operator/= (T val)
    {
      re() /= val; 
      im() /= val; 
    }


    /** \brief returns the sum of two complex numbers. 
     *
     * \return the sum of the complex numbers to the left and right of the operator. */
    template <class T> inline ComplexNumber<T> ComplexNumber<T>::operator+ (const ComplexNumber<T>& c) const 
    {
      ComplexNumber<T> ret (re()+c.re(), im()+c.im());
      return (ret); 
    }


    /** \brief returns the subtraction of two complex numbers. 
     *
     * \return the subtraction of the complex numbers to the left and right of the operator. */
    template <class T> inline ComplexNumber<T> ComplexNumber<T>::operator- (const ComplexNumber<T>& c) const
    {
      ComplexNumber<T> ret (re()-c.re(), im()-c.im());
      return (ret); 
    }


    /** \brief returns the product of two complex numbers. 
     *
     * \return the product of the complex numbers to the left and right of the operator. */
    template <class T> inline ComplexNumber<T> ComplexNumber<T>::operator* (const ComplexNumber<T>& c) const
    {
      ComplexNumber<T> ret (re()*c.re() - im()*c.im(), re()*c.im() + im()*c.re());
      return (ret); 
    }


    /** \brief adds a complex number. 
     *
     * adds the complex number on the right of the operator to the complex number on the left of the operator*/
    template <class T> inline void ComplexNumber<T>::operator+= (const ComplexNumber<T>& c) 
    {
      re() += c.re(); 
      im() += c.im(); 
    }


    /** \brief subtracts a complex number. 
     *
     * subtracts the complex number on the right of the operator from the complex number on the left of the operator*/
    template <class T> inline void ComplexNumber<T>::operator-= (const ComplexNumber<T>& c)
    {
      re() -= c.re(); 
      im() -= c.im(); 
    }


    /** \brief multiplies by a complex number. 
     *
     * multiplies the complex number on the left of the operator by the complex number on the right of the operator*/
    template <class T> inline void ComplexNumber<T>::operator*= (const ComplexNumber<T>& c)
    {
      T r = re()*c.re() - im()*c.im();
      im() = re()*c.im() + im()*c.re();
      re() = r;
    }

    /** \brief returns the complex conjugate of the number. 
     *
     * \return the complex conjugate of the number. */
    template <class T> inline ComplexNumber<T> ComplexNumber<T>::conj () const
    {
      ComplexNumber<T> ret (re(), -im()); 
      return (ret); 
    }


    /** \brief returns a normalised copy of the complex number. 
     *
     * \return a copy of the complex number, normalised to unit modulus. */
    template <class T> inline ComplexNumber<T> ComplexNumber<T>::norm () const
    {
      return ((*this)/mod());
    }


    /** \brief returns a pointer to the data. 
     *
     * used to access the underlying complex data directly.
     * \return a pointer to the data array holding the real and imaginary parts of the complex number. 
     * \warning only access the complex data directly if the existing functions are not adequate. */
    template <class T> inline T* ComplexNumber<T>::pointer ()
    {
      return (data); 
    }


    /** \brief print out complex number to the specified stream. 
     *
     * this operator prints out the complex number onto the specified C++ output stream. */
    template <class T> inline std::ostream& operator<< (std::ostream& stream, const ComplexNumber<T>& c) 
    {
      stream << c.re() << " + " << c.im() << "i"; 
      return (stream); 
    }

  }
}



#endif

