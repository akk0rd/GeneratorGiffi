#include<iostream>
#include<omp.h>
#include<cmath>
#include<bitset>
#include<fstream>
#include<string>
#include<vector>
#include<bitset>

const unsigned N1 = 265;
const unsigned C1 = 83;
const unsigned N2 = 272;
const unsigned C2 = 85;
const unsigned N = 2048;

std::string input(std::string inp);
std::pair<double, double> reduceNum();
std::pair<double, double> reduceNum2();

std::bitset<N> takeBit(std::string text);
size_t minV(size_t fir, size_t sec);
size_t statistics(std::bitset<N> end, std::bitset<N1> first);
std::vector<uint> generator(std::bitset<30> L, std::bitset<N> G);
std::vector<uint> generator(std::bitset<31> L, std::bitset<N> G);
void endFound(std::vector<uint> candidatsL1, std::vector<uint> candidatsL2,std::bitset<N> G,std::bitset<30> L1,std::bitset<31> L2,std::bitset<32> L3);

int main() {
    std::string inp = input("input.txt");

    std::pair<double,double> rN1,rN2;
    rN1=reduceNum();
    rN2=reduceNum2();
    std::cout<<"Number of digits for L1 "<<rN1.first<<" and threshold "<<rN1.second<<std::endl;
    std::cout<<"Number of digits for L2 "<<rN2.first<<" and threshold "<<rN2.second<<std::endl;
    std::bitset<N> vec = takeBit(inp);
    std::bitset<30> L1;
    L1[30] = 1, L1[6] = 1, L1[4] = 1, L1[1] = 1, L1[0] = 1;
    std::bitset<31> L2;
    L2[31] = 1, L2[3] = 1, L2[0] = 1;
    std::bitset<32> L3;
    L3[32] = 1, L3[7] = 1, L3[5] = 1, L3[3] = 1, L3[2] = 1, L3[1] = 1, L3[0] = 1;

    std::vector<uint> a;
    std::vector<uint> b;
    a = generator(L1, vec);
    b = generator(L2, vec);
    endFound(a,b,vec,L1,L2,L3);
    return EXIT_SUCCESS;
}

std::string input(std::string inp) {
    std::ifstream file;
    file.open(inp);
    std::string text;
    std::getline(file, text);
    return text;
}

std::pair<double,double> reduceNum() {
    float p1 = 0.25;
    float p2 = 0.5;
    float qalph = 2.326348;
    float qbetta=6.120756058;
    float k1 = qalph*sqrt(p1*(1 - p1));
    float k2 = qbetta*sqrt(p2*(1 - p2));
    double Nz = std::pow(((k1 + k2) / (p1 - p2)), 2);
    double C = Nz*p2 - sqrt(Nz*p2*(1 - p2))*qbetta;
    std::pair<double, double> rez(Nz, C);
    return rez;
}
std::pair<double,double> reduceNum2() {
    float p1 = 0.25;
    float p2 = 0.5;
    float qalph = 2.326348;
    float qbetta=6.230259914;
    float k1 = qalph*sqrt(p1*(1 - p1));
    float k2 = qbetta*sqrt(p2*(1 - p2));
    double Nz = std::pow(((k1 + k2) / (p1 - p2)), 2);
    double C = Nz*p2 - sqrt(Nz*p2*(1 - p2))*qbetta;
    std::pair<double, double> rez(Nz, C);
    return rez;
}

std::bitset<N> takeBit(std::string text) {
    std::bitset<N> vec;
    char a;
    for (size_t i = 0;i < N;i++) {
        char a = text[i];
        vec[N-1-i] = (bool)atoi(&a);
    }
    return vec;
}
std::vector<uint> generator(std::bitset<30> L, std::bitset<N> G) {
    std::ofstream f;
    f.open("end.txt");
    std::vector<uint> trueBits;
#pragma omp parallel for
    for (int i = 0;i < (int)pow(2, 30);i++) {
        int R=0;
        for(size_t k = 0; k< 30;k++){
            R+=L[k] ^ G[N-k];
        }
        std::bitset<30> tmp(i);
        for (size_t j = 0;j < N1;j++){
            std::bitset<30> cnt = tmp & L;
            tmp <<= 1;
            int p = cnt.count() % 2;
            tmp[0] = p;
            R+=p ^ G[N-j-30];
        }
        if(R<C1){
            trueBits.push_back(i);
        }
    }
    f.close();
    return trueBits;
}
std::vector<uint> generator(std::bitset<31> L, std::bitset<N> G) {
    std::ofstream f;
    f.open("endL2.txt");
    std::vector<uint> trueBits;
#pragma omp parallel for
    for (int i = 0;i < (int)pow(2, 31);i++) {
        int R=0;
        for(size_t k = 0; k< 31;k++){
            R+=L[k] ^ G[N-k];
        }
        std::bitset<31> tmp(i);
        for (size_t j = 0;j < N2;j++){
            std::bitset<31> cnt = tmp & L;
            tmp <<= 1;
            int p = cnt.count() % 2;
            tmp[0] = p;
            R+=p ^ G[N-j-31];
        }
        if(R<C2){
            trueBits.push_back(i);
        }
    }
    f.close();
    return trueBits;
}

void endFound(std::vector<uint> candidatsL1, std::vector<uint> candidatsL2,std::bitset<N> G,std::bitset<30> L1,std::bitset<31> L2,std::bitset<32> L3){

    bool Flag = false;
#pragma omp parallel for
    for (unsigned long int i = 0; i < pow(2,32); i++) {

        for (int k = 0; k < candidatsL1.size(); k++)
        {

            for (int l = 0; l < candidatsL2.size(); l++)
            {
                std::bitset<30> stateL1(candidatsL1[k]);
                std::bitset<31> stateL2(candidatsL2[l]);
                std::bitset<32> stateL3(i);

                bool answer = true;

                for (int i = 0; i < N; i++)
                {
                    if ((stateL3[32]) == 1) {
                        if ((stateL1[30]) != G[N-i]) {
                            answer = false;
                            break;
                        }

                    } else {
                        if (stateL2[31] != G[N-i]) {
                            answer = false;
                            break;
                        }

                    }
                    std::bitset<30> tmp1 = stateL1 & L1;
                    stateL1<<=1;
                    stateL1[0]=tmp1.count()%2;

                    std::bitset<31> tmp2 = stateL2 & L2;
                    stateL2<<= 1;
                    stateL2[0]= tmp2.count()%2;

                    std::bitset<32> tmp3 = stateL3 & L3;
                    stateL3<<=1;
                    stateL3[0] = tmp3.count()%2;
                }

                if (answer) {
                    std::cout << "L1: " << candidatsL1[k] << std::endl;
                    std::cout << "L2: " << candidatsL2[l] << std::endl;
                    std::cout << "L3: " << i << std::endl;
                    Flag = true;
                    break;
                }
            }

            if (Flag) break;
        }

        if (Flag) break;
    }
}