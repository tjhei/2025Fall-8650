// Naive recursive version of factorial()
//
// What happens if you turn optimizations to -O2?

#include <iostream>

int factorial(int num) {
    if (num <= 1) 
        return 1;

    int result = num*factorial(num-1);
    return result;
}

int main()
{
  // increase to 1000000 to see a crash:
  std::cout << factorial(10) << std::endl;
}

/*
   rsp: stack pointer (top)
   rbp: base pointer (bottom of current function)
   [rbp-X]: local variables
   [rbp+X]: function arguments

enter:
    push    rbp      - save caller's base pointer
    mov     rbp, rsp - copy stack pointer rsp into rbp
    sub     rsp, 16  - make 16 bytes of stack space

recurse:
    call             - put return address on stack, jump to function

exit:
    leave            - set rsp to rbp, pop top of stack, assign to rbp
    ret              - pop return address from stack, jump there


*/
