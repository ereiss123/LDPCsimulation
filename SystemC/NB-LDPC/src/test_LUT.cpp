#include "itpp/comm/galois.h"
#include <vector>
#include <iostream>

void print_addition_table(int q);
void print_multiplication_table(int q);

int main(void){
    int q = 8;
    std::cout << "Starting program" << std::endl;
    itpp::GF z = itpp::GF(q,-1); // return 0th element
    itpp::GF z2 = itpp::GF(q,2); // return 0th element
    std::cout << "Field set" << std::endl;
    std::cout << z.get_size() << std::endl;
    std::cout << z.get_value() << std::endl;
    itpp::GF sum = z+z2;
    std::cout << sum.get_size() << std::endl;
    std::cout << sum.get_value() << std::endl;
    // std::cout << (z+z2) << std::endl;
    std::vector<std::vector<std::pair<int,int> > >  LUT;
  // iterate over every symbol combination
  for(int i = -1; i < q-1 ; i++)
  {
    itpp::GF GF_i(q,i);
    std::vector<std::pair<int,int> > idxs;

    for(int j = -1; j < q-1; j++)
    {
      itpp::GF GF_j(q,j);
      if(GF_j + GF_i == z)
      {
        std::cout << "Testing " << GF_i.get_value() << ' ' << GF_j.get_value() << std::endl;
        std::cout << "Equal 0 element? "<< ((GF_j + GF_i == z)?"True":"False")  << std::endl;
        std::pair<int,int> elements(i,j);
        idxs.emplace_back(elements);
      }
    }
    LUT.emplace_back(idxs);
  }
    print_addition_table(4);
    std::cout << "\n" << std::endl;
    print_addition_table(8);
    std::cout << "\n" << std::endl;
    print_multiplication_table(4);
    std::cout << "\n" << std::endl;
    print_multiplication_table(8);

    return 0;
}

void generate_LUT(int q)
{
  itpp::GF z = itpp::GF(q); // return 0th element
  // iterate over every symbol combination
  for(int i = 0; i < q ; i++)
  {
    itpp::GF GF_i(q,i);
    std::vector<std::pair <int, int> > idxs;

    itpp::GF i_gf = itpp::GF(p.alist.q,i);
    for(int j = 0; j < p.alist.q; j++)
    {
      itpp::GF GF_j(p.alist.q,j);
      if(GF_j + GF_i == z)
      {
        std::pair<int,int> elements(i+1,j+1);
        idxs.emplace_back(elements);
      }
    }
    p.LUT.emplace_back(idxs);
  }
}

void print_addition_table(int q){
    std::cout << "+ |";
    // Print addition table
    for(int i = -1; i < q-1 ; i++)
    {
        if(i == -1)
        {
            std::cout << ' ' << i ;
        }
        else
        {
            std::cout << "  " << i ;
        }
    }

    std::cout << std::endl;
    for (size_t i = 0; i < (q*3)+3; i++)
    {
        std::cout << '-';
    }


  std::cout << std::endl;
  for(int i = -1; i < q-1 ; i++)
  {
    itpp::GF GF_i(q,i);
    std::vector<std::pair<int,int> > idxs;

    for(int j = -2; j < q-1; j++)
    {
        if(j == -2)
        {
            if(i == -1)
            {
                std::cout << i << '|' ;
            }
            else
            {
                std::cout <<  i << " |" ;
            }
        }
        else
        {
            itpp::GF GF_j(q,j);
            itpp::GF sum = GF_i + GF_j;
            if(sum.get_value() == -1)
            {
                std::cout << ' ' << sum.get_value();
            }
            else
            {
                std::cout << "  " << sum.get_value();
            }
        }
    }
    std::cout << std::endl;
  }
}
void print_multiplication_table(int q){
    std::cout << "* |";
    // Print addition table
    for(int i = -1; i < q-1 ; i++)
    {
        if(i == -1)
        {
            std::cout << ' ' << i ;
        }
        else
        {
            std::cout << "  " << i ;
        }
    }

    std::cout << std::endl;
    for (size_t i = 0; i < (q*3)+3; i++)
    {
        std::cout << '-';
    }


  std::cout << std::endl;
  for(int i = -1; i < q-1 ; i++)
  {
    itpp::GF GF_i(q,i);
    std::vector<std::pair<int,int> > idxs;

    for(int j = -2; j < q-1; j++)
    {
        if(j == -2)
        {
            if(i == -1)
            {
                std::cout << i << '|' ;
            }
            else
            {
                std::cout <<  i << " |" ;
            }
        }
        else
        {
            itpp::GF GF_j(q,j);
            itpp::GF prod = GF_i * GF_j;
            if(prod.get_value() == -1)
            {
                std::cout << ' ' << prod.get_value();
            }
            else
            {
                std::cout << "  " << prod.get_value();
            }
        }
    }
    std::cout << std::endl;
  }
}