#include<iostream>
class A {
   public: 
   int b;
   A() 
   { 
       b = 5; 
    }
};
int test() {
   A a = A();
   A* x = &a;
   std::cout << "a.b = " << a.b << "\n";
   std::cout << "x->b = " << x->b << "\n";
   return 0;
}