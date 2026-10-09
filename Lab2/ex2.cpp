#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int main()
{
    string S = "TAAGATTTCCTAGGT";
    int f[256] = {0};

    for(int i = 0; i < S.size(); i++)
    {
        f[S[i]]++;
    }

    double na = 0.05;
    double gcper = 100.0 * (f['C']+f['G']) / S.size();
    double Tm2 = 81.5 + 16.6 * (log10(na)) + 0.41 * gcper - 600.0 / S.size();

    cout<<"Actual Tm: "<<Tm2<<"\n";
    cout<<"Freq of A: "<<f['A']<<"\n";
    cout<<"Freq of T: "<<f['T']<<"\n";
    cout<<"Freq of C: "<<f['C']<<"\n";
    cout<<"Freq of G: "<<f['G']<<"\n";

    return 0;
}