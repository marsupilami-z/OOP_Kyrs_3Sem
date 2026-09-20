#pragma once

#include <string>
#include <functional>
#include <stdexcept>
#include <cmath>

template <typename T>
using MF = std::function<T(T)>;
using std::string;

template <typename T>
class MathFunction 
{
public:
  MathFunction(MF<T> f, MF<T> df, string str_expr);
  MathFunction(MF<T> f, string str_expr);

  T eval(T x) const;
  T deriv(T x) const;

  bool has_deriv() const;
  const string& str_expr() const;

private:
  MF<T> m_func;
  MF<T> m_dfunc;

  string m_str_expr;
};

template <typename T>
class RootFinder
{
public:
  RootFinder(const MathFunction<T>& func, T eps, int max_iter);

  virtual T solve(T a, T b) = 0;
  const string& method_name() const;
  int real_iterations() const;
  
  virtual ~RootFinder() = default;

protected:
  const MathFunction<T>& m_func; 
  T m_eps;
  int m_max_iter;
  int m_real_iter;
  string m_method_name;
};

template <typename T>
class SimpleIteration : public RootFinder<T>
{
public:
  SimpleIteration(const MathFunction<T>& func, MF<T> phi, T eps, int max_iter);

  T solve(T a, T b) override;

private:
  MF<T> m_phi;
};

template <typename T>
class Newton : public RootFinder<T>
{
public:
  Newton(const MathFunction<T>& func, T eps, int max_iter);
  T solve(T a, T) override;
};

template <typename T>
class Secant : public RootFinder<T>
{
public:
  Secant(const MathFunction<T>& func, T eps, int max_iter);
  T solve(T a, T) override;
};

template <typename T>
class Bisection : public RootFinder<T>
{
public:
  Bisection(const MathFunction<T>& func, T eps, int max_iter);
  T solve(T a, T) override;  
};

// --- ///

template <typename T>
MathFunction<T>::MathFunction(MF<T> f, MF<T> df, string str_expr) 
  : m_func(f), m_dfunc(df), m_str_expr(str_expr)
{

}

template <typename T>
MathFunction<T>::MathFunction(MF<T> f, string str_expr) 
  : m_func(f), m_dfunc(nullptr), m_str_expr(str_expr) 
{
  
}

template <typename T>
T MathFunction<T>::eval(T x) const { return m_func(x); }

template <typename T>
T MathFunction<T>::deriv(T x) const 
{ 
  if(!m_dfunc) throw std::runtime_error("dfunc not provided");
  return m_dfunc(x); 
}

template <typename T>
bool MathFunction<T>::has_deriv() const 
{
  if(m_dfunc) return true;
  return false;
}

template <typename T>
const string& MathFunction<T>::str_expr() const { return m_str_expr; }

template <typename T>
RootFinder<T>::RootFinder(const MathFunction<T>& func, T eps, int max_iter)
  : m_func(func), m_eps(eps), m_max_iter(max_iter), m_real_iter(0)
{
  // --- Недопустимо
  if(eps <= 0) 
    throw std::invalid_argument("eps should be > 0");
  if(max_iter <= 0)
    throw std::invalid_argument("max iter should be > 0");
  // ---
}

template <typename T>
int RootFinder<T>::real_iterations() const { return m_real_iter; } 

template <typename T>
const string& RootFinder<T>::method_name() const { return m_method_name; }

template <typename T>
SimpleIteration<T>::SimpleIteration(const MathFunction<T>& func, MF<T> phi, T eps, int max_iter)
  : RootFinder<T>(func, eps, max_iter), m_phi(phi)
{
  this->m_method_name = "Simple Iteration";

  if(!phi) throw std::invalid_argument("phi should exist"); 
}

template <typename T>
T SimpleIteration<T>::solve(T a, T)
{
  T x0{ a };
  T h{ 1e-7 };

  T phi_1{ m_phi(x0+h) - m_phi(x0-h) };
  T dphi{ phi_1 / (2*h) };

  if(std::abs(dphi) > 1) throw std::runtime_error("dphi > 1 error");

  this->m_real_iter = 0;
  T x_prev{a};
  T x_next;

  while( this->m_real_iter < this->m_max_iter )
  {
    x_next = m_phi(x_prev);
    this->m_real_iter++;

    if( std::abs(x_next - x_prev) < this->m_eps ) { return x_next; }
    else { x_prev = x_next; }
  }

  throw std::runtime_error("SimpleIteration: did not converge");
}

template <typename T>
Newton<T>::Newton(const MathFunction<T>& func, T eps, int max_iter)
  : RootFinder<T>(func, eps, max_iter)
{
  this->m_method_name = "Newton";
}

template <typename T>
T Newton<T>::solve(T a, T)
{
  this->m_real_iter = 0;
  T x_prev{a};
  T x_next;

  T h{1e-7};

  T df;

  bool has_deriv{ this->m_func.has_deriv() };

  while(this->m_real_iter < this->m_max_iter)
  {
    if(has_deriv)
    {
      df = this->m_func.deriv(x_prev);
    }
    else
    {
      df = (this->m_func.eval(x_prev+h) - this->m_func.eval(x_prev-h)) / (2*h);
    }

    if( std::abs(df) < 1e-15 )
    {
      throw std::runtime_error("Newton: derivative is near zero at x");
    }

    x_next = x_prev - this->m_func.eval(x_prev) / df;

    this->m_real_iter++;

    if( std::abs(x_next - x_prev) < this->m_eps) { return x_next; }
    x_prev = x_next;
  }

  throw std::runtime_error("Newton: did not converge");
}

template <typename T>
Secant<T>::Secant(const MathFunction<T>& func, T eps, int max_iter)
  : RootFinder<T>(func, eps, max_iter)
{
  this->m_method_name = "Secant";
}

template <typename T>
T Secant<T>::solve(T a, T b)
{
  if(a == b) throw std::invalid_argument("Secant: a = b");

  this->m_real_iter = 0;

  T f_0{this->m_func.eval(a)};
  T f_1{this->m_func.eval(b)};

  while(this->m_real_iter < this->m_max_iter)
  {
    T denom{f_1 - f_0};

    if(std::abs(denom) < 1e-15) throw std::runtime_error("Secant: denom zero");

    T x_next{ b - (f_1 * (b-a)) / denom };
    this->m_real_iter++;

    if(std::abs(x_next - b) < this->m_eps) return x_next;

    a = b;
    f_0 = f_1;
    b = x_next;
    f_1 = this->m_func.eval(b);
  }

  throw std::runtime_error("Secant: no converge");
}

template <typename T>
Bisection<T>::Bisection(const MathFunction<T>& func, T eps, int max_iter)
  : RootFinder<T>(func, eps, max_iter)
{
  this->m_method_name = "Bisection";
}

template <typename T>
T Bisection<T>::solve(T a, T b)
{
  if(a >= b)
    throw std::invalid_argument("Bisection: a should be less than b");
    
  T f_a{ this->m_func.eval(a) };
  T f_b{ this->m_func.eval(b) };

  if (std::abs(f_a) < this->m_eps) return a;
  if (std::abs(f_b) < this->m_eps) return b;

  if (f_a * f_b > 0)
    throw std::invalid_argument("Bisection: no root here");

  this->m_real_iter = 0;

  while(this->m_real_iter < this->m_max_iter)
  {
    T c{ (a + b) / 2 };
    T f_c{ this->m_func.eval(c) };

    this->m_real_iter++;

    if((b-a) / 2 < this->m_eps) return c;

    if(f_a*f_c < 0)
    {
      b = c;
    }
    else
    {
      a = c;
      f_a = f_c;
    }
  }

  throw std::runtime_error("Bisection: no converge");
}