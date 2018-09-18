/*
 * (C) 2018 University of Bristol. See License.txt
 *
 * Demonstrate external client inputing and receiving outputs from a SPDZ process, 
 * following the protocol described in https://eprint.iacr.org/2015/1006.pdf.
 *
 * Provides a client to bankers_bonus.mpc program to calculate which banker pays for lunch based on
 * the private value annual bonus. Up to 8 clients can connect to the SPDZ engines running 
 * the bankers_bonus.mpc program.
 *
 * Each connecting client:
 * - sends a unique id to identify the client
 * - sends an integer input (bonus value to compare)
 * - sends an integer (0 meaining more players will join this round or 1 meaning stop the round and calc the result).
 *
 * The result is returned authenticated with a share of a random value:
 * - share of winning unique id [y]
 * - share of random value [r]
 * - share of winning unique id * random value [w]
 *   winning unique id is valid if ∑ [y] * ∑ [r] = ∑ [w]
 * 
 * No communications security is used. 
 *
 * To run with 2 parties / SPDZ engines:
 *   ./Scripts/setup-online.sh to create triple shares for each party (spdz engine).
 *   ./compile.py bankers_bonus
 *   ./Scripts/run-online bankers_bonus to run the engines.
 *
 *   ./bankers-bonus-client.x 123 2 100 0
 *   ./bankers-bonus-client.x 456 2 200 0
 *   ./bankers-bonus-client.x 789 2 50 1
 *
 *   Expect winner to be second client with id 456.
 */

#include "Math/gfp.h"
#include "Math/gf2n.h"
#include "Networking/sockets.h"
#include "Tools/int.h"
#include "Math/Setup.h"
#include "Auth/fake-stuff.h"

#include <sodium.h>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <cstdlib>



const int NUM_ROWS = 10;
const int NUM_COLUMNS = 10;

// Send the private inputs masked with a random value.
// Receive shares of a preprocessed triple from each SPDZ engine, combine and check the triples are valid.
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


// Generates matrix with a finish
vector<gfp> generateMatrix(int finish) {
    cout << "Generating Matrix" << endl;
    vector<gfp> matrix(1 + NUM_ROWS * NUM_COLUMNS);
    gfp finish_2 = finish;
    matrix[0] = finish_2;
    for (int i = 1; i < NUM_ROWS * NUM_COLUMNS + 1; i++) {
        gfp val = rand();
        matrix[i] = val;
    }
    for(unsigned int i=0; i<matrix.size(); ++i)
        cout << matrix[i] << ',';

    return matrix;
}


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


vector<gfp> receive_result(vector<int>& sockets, int nparties)
{
    cout << "Receiving matrix" << endl;
    vector<gfp> output_values(NUM_ROWS * NUM_COLUMNS);    
    octetStream os;
    for (int i = 0; i < NUM_ROWS * NUM_COLUMNS; i++)
    {
        output_values[i] = receive_one_result(sockets, nparties);
    }

    return output_values;
}





int main(int argc, char** argv)
{
    srand(time(NULL));
    clock_t start;
    double duration;
    int port_base = 14000;
    int nparties = 2;
    int finish;
    string host_name = "localhost";

    if (argc < 1) {
        cout << "Please provide client id" << endl;
        exit(0);
    }


    finish = atoi(argv[0]);
    //nparties = atoi(argv[1]);

    // init static gfp
    string prep_data_prefix = get_prep_dir(nparties, 128, 40);
    initialise_fields(prep_data_prefix);


    vector<int> sockets(nparties);
    for (int i = 0; i < nparties; i++)
    {
        set_up_client_socket(sockets[i], host_name.c_str(), port_base + i);
    }
    cout << "Finish setup socket connections to SPDZ engines." << endl;
    vector<gfp> matrix = generateMatrix(finish);


    start = clock();
    // Run the commputation
    send_private_inputs(matrix, sockets, nparties);
    cout << "Sent private inputs to each SPDZ engine, waiting for result..." << endl;

    // Get the result back (client_id of winning client)
    vector<gfp> result = receive_result(sockets, nparties);

    duration = (clock() - start ) / (double) CLOCKS_PER_SEC;
    printf(" Took %f seconds for matrix multiplication", duration);
    for (int i = 0; i < NUM_ROWS * NUM_COLUMNS; i++) {
        cout << result[i] << " , ";
    }
    for (int i = 0; i < nparties; i++) {
        close_client_socket(sockets[i]);
    }
    return 0;
}
