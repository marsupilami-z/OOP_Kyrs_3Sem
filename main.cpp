/* задача: Поиск корней непрерывной функции
 *
 * Уравнение вида f(x) = 0
 *
 * Реализовать следующие численные методы уточнения корней:
 * 1) Метод последовательных приближений (метод итераций)
 * 2) Метод Ньютона (метод касательных)
 * 3) Метод секущих (метод хорд)
 * 4) Метод половинного деления (метод дихотомии)
 *
 * Задать несколько функций для тестирования.
 * Реализовать меню выбора метода
 */

/*
 * В файле main.cpp будут исключительно тесты и интерфейс.
 * В файлах root_find.h реализация алгоритмов
 */

#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>

#include "root_find.h"

int main()
{
  unsigned int user_func_choose;
  unsigned int user_method_choose;

  while (true)
  {
    std::cout << "Choose function to solve: \n"
              << "1. x^2 - 2 \n"
              << "2. cos(x) - x \n"
              << "3. exp(x) - 3x \n"
              << "4. x^3 - x - 2 \n\n";

    std::cin >> user_func_choose;

    if (user_func_choose >= 1 && user_func_choose <= 4) break;

    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  MathFunction<double> func{nullptr, nullptr, ""};
  MF<double> phi = nullptr;

  switch (user_func_choose)
  {
  case 1:
    func = MathFunction<double>(
        [](double x) { return x*x - 2.0; },
        [](double x) { return 2.0 * x; },
        "x^2 - 2"
    );
    phi = [](double x) { return (x + 2.0 / x) / 2.0; };
    break;

  case 2:
    func = MathFunction<double>(
        [](double x) { return std::cos(x) - x; },
        [](double x) { return -std::sin(x) - 1.0; },
        "cos(x) - x"
    );
    phi = [](double x) { return std::cos(x); };
    break;

  case 3:
    func = MathFunction<double>(
        [](double x) { return std::exp(x) - 3.0 * x; },
        [](double x) { return std::exp(x) - 3.0; },
        "exp(x) - 3x"
    );
    phi = [](double x) { return std::exp(x) / 3.0; };
    break;

  case 4:
    func = MathFunction<double>(
        [](double x) { return x*x*x - x - 2.0; },
        [](double x) { return 3.0 * x*x - 1.0; },
        "x^3 - x - 2"
    );
    phi = [](double x) { return std::cbrt(x + 2.0); };
    break;
  }

  while (true)
  {
    std::cout << "Choose method: \n"
              << "1. Simple Iterations \n"
              << "2. Newton \n"
              << "3. Secant \n"
              << "4. Bisection \n\n";

    std::cin >> user_method_choose;

    if (user_method_choose >= 1 && user_method_choose <= 4) break;

    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  double a, b, eps;
  int max_iter;

  std::cout << "Enter a: ";         
  std::cin >> a;

  std::cout << "Enter b: ";         
  std::cin >> b;

  std::cout << "Enter eps: ";       
  std::cin >> eps;

  std::cout << "Enter max_iter: ";  
  std::cin >> max_iter;

  auto print_result = [&func](RootFinder<double>& method, double a, double b)
  {
    double root = method.solve(a, b);

    std::cout << "\nFunction: " << func.str_expr() << "\n"
              << "Method:   " << method.method_name() << "\n"
              << std::setprecision(12) << "Root:     " << root << "\n"
              << "Iterations: " << method.real_iterations() << "\n"
              << "|f(root)| = " << std::abs(func.eval(root)) << "\n";
  };


  switch (user_method_choose)
  {
    case 1:
    {
      SimpleIteration<double> method(func, phi, eps, max_iter);
      print_result(method, a, b);
      break;
    }

    case 2:
    {
      Newton<double> method(func, eps, max_iter);
      print_result(method, a, b);
      break;
    }

    case 3:
    {
      Secant<double> method(func, eps, max_iter);
      print_result(method, a, b);
      break;
    }

    case 4:
    {
      Bisection<double> method(func, eps, max_iter);
      print_result(method, a, b);
      break;
    }
  }

  return 0;
}