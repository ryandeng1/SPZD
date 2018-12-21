#include "Math/gfp.h"
#include "Math/gf2n.h"
#include "Networking/sockets.h"
#include "Tools/int.h"
#include "Math/Setup.h"
#include "Auth/fake-stuff.h"
#include "Eigen/Dense"
#include "json/json.hpp"
#include<gmp.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>

#include <sodium.h>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <cstdlib>

using namespace NTL;

// Add the private input value to triple[0] and send to each spdz engine.
void send_private_inputs(vector<gfp>& values, vector<int>& sockets, int nparties)
{
    int num_inputs = values.size();
    octetStream os;
    vector< vector<gfp> > triples(num_inputs, vector<gfp>(3));
    vector<gfp> triple_shares(3);

    // Receive num_inputs triples from SPDZ
    for (int j = 0; j < nparties; j++)
    {
        os.reset_write_head();
        os.Receive(sockets[j]);

        for (int j = 0; j < num_inputs; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                triple_shares[k].unpack(os);
                triples[j][k] += triple_shares[k];
            }
        }
    }

    // Check triple relations (is a party cheating?)
    for (int i = 0; i < num_inputs; i++)
    {
        if (triples[i][0] * triples[i][1] != triples[i][2])
        {
            cout << triples[i][0] << " , " << triples[i][1] << " , " << triples[i][2] << endl;
            cerr << "Incorrect triple at " << i << ", aborting\n";
            exit(1);
        }
    }

    os.reset_write_head();
    // Send inputs + triple[0], so SPDZ can compute shares of each value
    for (int i = 0; i < num_inputs; i++)
    {
        gfp y = values[i] + triples[i][0];
        y.pack(os);
        cout << y << " , ";
    }
    for (int j = 0; j < nparties; j++)
        os.Send(sockets[j]);
}

// Assumes that Scripts/setup-online.sh has been run to compute prime
void initialise_fields(const string& dir_prefix)
{
  int lg2;
  bigint p;

  string filename = dir_prefix + "Params-Data";
  cout << "loading params from: " << filename << endl;

  ifstream inpf(filename.c_str());
  if (inpf.fail()) { throw file_error(filename.c_str()); }
  inpf >> p;
  inpf >> lg2;
  inpf.close();

  gfp::init_field(p);
  gf2n::init_field(lg2);
}

// Based off of 
vector<gfp> keygen(int primeLength) {
    ZZ p, q;
    GenPrimePair(p, q, primeLength);
    ZZ n = p * q;
    ZZ g = n + 1;
    ZZ phi = (p - 1) * (q - 1);
    // LCM(p, q) = p * q / GCD(p, q);
   	// ZZ lambda = phi / GCD(p - 1, q - 1);
   	ZZ lambda = phi;
    ZZ mu = InvMod(lambda, modulus);

    vector<gfp> keys(4);
    keys[0] = n;
    keys[1] = g;
    keys[2] = lambda;
    keys[3] = mu;
    return keys;
}


/*
gfp receive_one_result(vector<int>& sockets, int nparties)
{
    vector<gfp> output_values(3);
    octetStream os;
    for (int i = 0; i < nparties; i++)
    {
        os.reset_write_head();
        os.Receive(sockets[i]);
        for (unsigned int j = 0; j < 3; j++)
        {
            gfp value;
            value.unpack(os);
            output_values[j] += value;            
        }
    }

    if (output_values[0] * output_values[1] != output_values[2])
    {
        cerr << "Unable to authenticate output value as correct, aborting." << endl;
        exit(1);
    }
    return output_values[0];
}



vector<double> receive_result(vector<int>& sockets, int nparties, int NUM_COLUMNS)
{
    cout << "Receiving matrix" << endl;
    vector<double> output_values(NUM_COLUMNS);    
    octetStream os;
    for (int i = 0; i < NUM_COLUMNS; i++)
    {
        //gfp gfp_val;
        gfp gfp_val = receive_one_result(sockets, nparties);


        const gfp gfp_regular = gfp_val;
        bigint val;
        to_bigint(val, gfp_regular);

        bigint val_negate;
        gfp_val.negate();
        to_bigint(val_negate, gfp_val);


        double converted_double = mpz_get_d(val.get_mpz_t()) / pow(2, 20);
        double converted_double_negate = -1 * mpz_get_d(val_negate.get_mpz_t()) / pow(2, 20);
        cout << "Converted double " << converted_double << endl;
        cout << "Converted double negative " << converted_double_negate << endl;
        if (abs(converted_double) < 10) {
            output_values[i] = converted_double;
        } else {
            output_values[i] = converted_double_negate;
        }
        //cout << "received " << output_values[i] << endl;
    }

    return output_values;
}
*/






void func(int argc, char** argv) {
    cout << argc;
    cout << argv;
    return;
}

int main(int argc, char** argv) {
    
    //srand(time(NULL));
    clock_t start;
    double duration;
    int port_base = 14000;
    //int nparties = 2;

    //Shift all numbers over by 20 bits
    int numShift = 20;

    //string host_names[] = {"ec2-52-39-162-238.us-west-2.compute.amazonaws.com", "ec2-34-223-215-198.us-west-2.compute.amazonaws.com", "ec2-23-20-124-131.compute-1.amazonaws.com", "ec2-52-73-142-253.compute-1.amazonaws.com"};

	string host_names[] = {"localhost", "localhost", "localhost", "localhost"};
    if (argc < 1) {
        cout << "Please provide client id" << endl;
        exit(0);
    }


    int nparties = atoi(argv[2]);
    // Default prime length to 128 bits
    int prime_length = 128;
    // Init
    string prep_data_prefix = get_prep_dir(nparties, 128, 128);
    initialise_fields(prep_data_prefix);




    vector<int> sockets(nparties);
    
    for (int i = 0; i < nparties; i++)
    {
        set_up_client_socket(sockets[i], host_names[i].c_str(), port_base + i);
    }
    cout << "Finish setup socket connections to SPDZ engines." << endl;
    
    start = clock();

    vector<gfp> keys = keygen(prime_length);
    send_private_inputs(values, sockets, nparties);
    cout << "Sent private inputs to each SPDZ engine, waiting for result..." << endl;

    // Get the result back (client_id of winning client)
    //vector<double> result = receive_result(sockets, nparties, cols);


    duration = (clock() - start ) / (double) CLOCKS_PER_SEC;
    printf(" Took %f seconds for Paillier Keygen", duration);
    
    for (int i = 0; i < nparties; i++) {
        close_client_socket(sockets[i]);
    }
    
    func(argc, argv);
}
