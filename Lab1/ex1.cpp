#include <iostream>
#include <string>

using namespace std;

int main()

{
string S = "TAAGATTTCCTAGGT";
int f[256] = {0};

for(int i = 0; i < S.size(); i++)
{
f[S[i]]++;
}

cout<<"Freq of A: "<<f['A']<<"\n";
cout<<"Freq of T: "<<f['T']<<"\n";
cout<<"Freq of C: "<<f['C']<<"\n";
cout<<"Freq of G: "<<f['G']<<"\n";

    return 0;
}