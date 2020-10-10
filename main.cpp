//
//  main.cpp
//  
//
//  Created by david keating on 6/13/19.
//

#include <stdio.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <climits>
#include <sstream>
#include <cmath>
#include <stack>
#include <random>
#include <chrono>

// critical Beta * J = log(1+sqrt(2)) / 2 for square lattice
// critical Beta * J = log(2+sqrt(2)) / 2 for hexagonal and triangular lattices

void ConfigToTxt(std::vector<int> &c, int M, int N, std::string filename) {
    std::ofstream outputFile(filename.c_str());
    
    
    for (int i = 0; i < N+4; i++) {
        for (int j = 0; j < M+4; j++) {
            outputFile<<c[i*(M+4)+j]<<" ";
            if (j==M+3) {
                outputFile<<"\n";
            }
        }
    }
    outputFile.close();
}

void Kawasaki(std::vector<int> &c, int M, int N, int Nsteps) {
    // Vertical edges: M+3 x N+4 but dont want to select the top and bottom two rows or the left or right two cols, so M-1 x N
    // Horizontal edges: M+4 x N+3 but dont want to select the left and right two cols or the top and bottom two rows, so M x N-1
    // Here we assume M,N >= 2
    std::random_device seed;
    std::default_random_engine generator(seed());
    std::uniform_int_distribution<int> randV(0,(M-1)*N-1);
    std::uniform_int_distribution<int> randH(0,M*(N-1)-1);
    std::uniform_real_distribution<> dist(0,1);
    
    int k=0;
    double BJ = .55;
    while (k < Nsteps) {
        float r = dist(generator);
        //std::cout<<r<<" "<<float((M-1)*N) / float(M*(N-1)+(M-1)*N) << std::endl;
        if (r < float((M-1)*N) / float(M*(N-1)+(M-1)*N)) {
            // Vertical
            int r1 = randV(generator); // = (M-1)*y+x
            int x = (r1 % (M-1)) + 2;
            int y = (r1-x+2)/(M-1) + 2;
            if (y==1) {
                std::cout<<"Tried ("<<x<<","<<y<<") with r="<<r1<<"."<<std::endl;
            }
            float r2 = dist(generator);
            if (c[(M+4)*y+x]*c[(M+4)*y+(x+1)] == -1) {
                // Swap
                // neighbors of (x,y): (x-1,y), (x,y-1), (x,y+1)
                // neighbors of (x+1,y): (x+2,y), (x+1,y-1), (x+1,y+1)
                // H0 = B*J [v(x,y)*(v(x-1,y) + v(x,y-1) + v(x,y+1)) + v(x+1,y)*(v(x+2,y) + v(x+1,y-1) + v(x+1,y+1))
                // H1 = B*J [v(x+1,y)*(v(x-1,y) + v(x,y-1) + v(x,y+1)) + v(x,y)*(v(x+2,y) + v(x+1,y-1) + v(x+1,y+1))
                // H1 - H0 = B*J [v(x,y)*(v(x+2,y) + v(x+1,y-1) + v(x+1,y+1) - v(x-1,y) - v(x,y-1) - v(x,y+1)) + v(x+1,y)*(v(x-1,y) + v(x,y-1) + v(x,y+1) - v(x+2,y) - v(x+1,y-1) - v(x+1,y+1))
                // if r2 < new/old = e^(H1 - H0) flip
             
                int DeltaH = (c[(M+4)*y+x] - c[(M+4)*y+(x+1)])*(c[(M+4)*y+(x+2)] + c[(M+4)*(y-1)+(x+1)] + c[(M+4)*(y+1)+(x+1)]) + (c[(M+4)*y+(x+1)]-c[(M+4)*y+x])*(c[(M+4)*y+(x-1)] + c[(M+4)*(y-1)+x] + c[(M+4)*(y+1)+x]);
                
                //std::cout<<"Try flip: "<<r2<<" "<<H1<<" "<<H0<<" "<<DeltaH<<std::endl;
                if (r2 < exp(BJ*DeltaH) ) {
                    c[(M+4)*y+x] *= -1;
                    c[(M+4)*y+(x+1)] *= -1;
                }
            }
        } else {
            // Horizontal
            int r1 = randH(generator); // = M*y+x
            int x = (r1 % M) + 2;
            int y = (r1-x+2)/M + 2;
            float r2 = dist(generator);
            if (c[(M+4)*y+x]*c[(M+4)*(y+1)+x] == -1) {
                // Swap
                // neighbors of (x,y): (x,y-1), (x-1,y), (x+1,y)
                // neighbors of (x,y+1): (x,y+2), (x-1,y+1), (x+1,y+1)
                int DeltaH = (c[(M+4)*y+x] - c[(M+4)*(y+1)+x])*(c[(M+4)*(y+2)+x] + c[(M+4)*(y+1)+x-1] + c[(M+4)*(y+1)+x+1]) + (c[(M+4)*(y+1)+x] - c[(M+4)*y+x])*(c[(M+4)*(y-1)+x] + c[(M+4)*y+x-1] + c[(M+4)*y+x+1]);
                
                //std::cout<<"Try flip: "<<r2<<" "<<DeltaH<<std::endl;
                if (r2 < exp(BJ*DeltaH)) {
                    c[(M+4)*y+x] *= -1;
                    c[(M+4)*(y+1)+x] *= -1;
                }
            }
        }
        k++;
    }
}

void IsingBasicEx(int M, int N, double ratio, int Nsteps, int iter) {
    // Runs Kawasaki dynamics on the Ising model on an NxM rectangle. The number of minus flips is fixed such at (num minus) = ratio * (N*M). Each iteration of the kawasaki dynamics attempts Nsteps number of flips. There are iter number of iterations. After each iteration a text file containing the spin configuration is saved to the folder ./output
    
    std::vector<int> vertices((M+4)*(N+4),1); // (M+4)*(N+4) rectangular domain with a + spin at every vertex. Two layers of + at every boundary.
    
    int spinFinal = floor(ratio*(N*M)); // How many - spins we want
    
    // randomly flip Floor[ratio*totalSpins] to -
    std::random_device seed;
    std::default_random_engine generator(seed());
    std::uniform_int_distribution<int> randM(0,M-1);
    std::uniform_int_distribution<int> randN(0,N-1);
    
    int spinCount = 0;
    while(spinCount < spinFinal) {
        int x = randM(generator);
        int y = randN(generator);
        if(vertices[(y+2)*(M+4)+(x+2)] != -1) {
            vertices[(y+2)*(M+4)+(x+2)] = -1;
            spinCount++;
            //std::cout<<x<<" "<<y<<std::endl;
        }
    }
//    for(int i=0; i<; ++i) {
//        for(int j=0; j<; ++j) {
//
//        }
//    }
    
    ConfigToTxt(vertices, M, N, "./output/test0.txt");
    
    int k = 0;
    while(k < iter) {
        Kawasaki(vertices, M, N, Nsteps);
        k++;
        std::cout<<"Completed "<<k<<" iterations."<<std::endl;
        std::string outputk = "./test";
        
        ConfigToTxt(vertices, M, N, "./output/test" + std::to_string(k) + ".txt");
    }
    
    ConfigToTxt(vertices, M, N, "./output/end_test.txt");
}



int main() {
    IsingBasicEx(100,100,0.8,1000000,1000);
    return 0;
}

