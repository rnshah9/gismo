/** @file mpi_example.cpp

    @brief Testing MPI with G+Smo

    Execute (eg. with 10 processes):
       mpirun -np 10 ./bin/mpi_example

    or provide a hosts file on a cluster:
       mpirun -hostfile hosts.txt ./bin/mpi_example

    If your cluster is using srun:
       srun -N 10 ./bin/mpi_example

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Hofer, R. Schneckenleitner
*/

#include <gismo.h>
#include <time.h>
#include <unistd.h>

double                    next_random

  ( void )

{
  static int  initialized = 0;
  int         next;

  if ( ! initialized )
  {
    int  my_rank;
    int  flag;

    MPI_Initialized ( &flag );

    if ( flag )
    {
      MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );

      srand ( (unsigned int) my_rank );
    }

    initialized = 1;
  }

  next = rand ();

  return ((double) next / (double) RAND_MAX);
}

using namespace gismo;

int main(int argc, char **argv)
{
  const int  N    = 100;
  const int  root = 0;
  const int  tag  = 1;

  double     number;
  double     max_num;
  double     recv_num;
  int        proc_count;
  int        my_rank;
  int        iiter;

  gsCmdLine cmd("An example for testing MPI with G+Smo.\n"
      "Execute (eg. with 10 processes):                                      "
      "  *  mpirun -np 10 ./bin/mpi_example\n"
      "or provide a hosts file on a cluster:                                 "
      "  *  mpirun -hostfile hosts.txt ./bin/mpi_example\n"
      "If your cluster is using srun:                                        "
      "  *  srun -N 10 ./bin/mpi_example"
  );
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  // Conditional compilation
#ifdef GISMO_WITH_MPI
  gsInfo << "Gismo was compiled with MPI support.\n";
#else
  gsInfo << "Gismo was compiled without MPI support.\n";
#endif


  // Initialize the MPI environment
  const gsMpi & mpi = gsMpi::init(argc, argv);

  // Get current wall time
  double wtime = mpi.wallTime();

  // Get the world communicator
  gsMpiComm comm = mpi.worldComm();
  MPI_Request req;

  //Get size and rank of the processor
  proc_count = comm.size();
  my_rank = comm.rank();

  GISMO_ASSERT(proc_count > 1,"At least two processes are required.\n");

  typedef std::tuple<int,int,double,double,bool,bool> tuple_t;

  std::queue<tuple_t> m_queue;
  std::queue<int> m_workers;
  double tstart = 0.0;
  double tend   = 1.0;
  int N0        = 10;
  double dt0    = (tend - tstart) / N0;
  gsMpiStatus status;

  int globalID = 0;



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // INIT
  ////////////////////////////////////////////////////////////////////////////////////////////////
  std::map<double,double> m_map;

  int njobs = 0;

  double x,y;
  bool message = true;

  tuple_t send, receive;

  if ( my_rank == 0 )
  {
    gsInfo<<"Making queue ...\n";
    gsInfo<<"level\tid\tvalue\n";
    for (int l = 0; l!= 1; l++)
    {
      for (int k = 0; k!= std::pow(2,l)*N0+1; k++)
      {
        //            lvl,ID,x,y,refine,stop
        tuple_t tuple{l,k,std::pow(2,-l)*dt0*k,NULL,0,false};
        // tuple.level = l;
        // tuple.ID = k;
        // tuple.x = std::pow(2,-l)*dt0*k;
        // tuple.refine = 0;
        // tuple.stop = false;
        m_queue.push(tuple);
        globalID++;
        gsInfo<<l<<"\t"<<k<<"\t"<<std::pow(2,-l)*dt0*k<<"\n";
      }
    }
    gsInfo<<"finished.\n";

    printf ("Adding workers ...\n");
    for (int w = 1; w!=proc_count; w++)
      m_workers.push(w);

    while (!m_queue.empty() && !m_workers.empty())
    {
      // push
      send = m_queue.front();
      m_queue.pop();
      // Question 1a: How to send the tuple <int,int,double> to slave?
      // Question 2 : How to send the job to ANY slave (i.e. the first available one)

      gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to    "<<m_workers.front()<<": ID: "<<std::get<1>(send)<<"; x = "<<std::get<2>(send)<<"; y = "<<std::get<3>(send)<<"; level = "<<std::get<0>(send)<<"\n";
      comm.isend(&send, 1, m_workers.front(),&req,tag);
      m_workers.pop();
      njobs++;
    }

    while (njobs > 0)
    {
      gsInfo<<njobs<<" job(s) running\n";
      comm.recv(&receive,1,MPI_ANY_SOURCE,tag,&status);
      gsInfo<<"[MPI process "<<my_rank<<"] Received a job from "<<status.MPI_SOURCE<<": ID: "<<std::get<1>(receive)<<"; x = "<<std::get<2>(receive)<<"; y = "<<std::get<3>(receive)<<"; level = "<<std::get<0>(receive)<<"\n";
      njobs--;
      m_workers.push(status.MPI_SOURCE);

      if (y < 0.75 && y > 0.25)
      {
        gsInfo<<"[MPI process "<<my_rank<<"] special!"<<"\n";
        // m_queue.push({l,k,std::pow(2,-l)*dt0*k})
      }

      m_map[std::get<2>(receive)] = std::get<3>(receive);

      while (!m_queue.empty() && !m_workers.empty())
      {
        // push
        send = m_queue.front();
        m_queue.pop();
        // Question 1a: How to send the tuple <int,int,double> to slave?
        // Question 2 : How to send the job to ANY slave (i.e. the first available one)

        gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to    "<<m_workers.front()<<": ID: "<<std::get<1>(send)<<"; x = "<<std::get<2>(send)<<"; y = "<<std::get<3>(send)<<"; level = "<<std::get<0>(send)<<"\n";
        comm.isend(&send, 1, m_workers.front(),&req,tag);
        m_workers.pop();
        njobs++;
      }
    }

    std::get<5>(send) = true;
    for (int w = 1; w!=proc_count; w++)
      comm.isend(&send, 1,w,&req,tag);
  }
  else
  {
    while (true)
    {
      tuple_t tuple;
      comm.recv(&tuple,1,0,tag,MPI_STATUS_IGNORE);

      if (std::get<5>(tuple))
      {
        gsInfo<<"[MPI process "<<my_rank<<"] I have to stop!!"<<"\n";
        break;
      }
      gsInfo<<"[MPI process "<<my_rank<<"] Received a job from "<<0<<": ID: "<<std::get<1>(tuple)<<"; x = "<<std::get<2>(tuple)<<"; y = "<<std::get<3>(tuple)<<"; level = "<<std::get<0>(tuple)<<"\n";


      std::get<3>(tuple) = std::pow(std::get<2>(tuple),2);
      // Question 1b: How to receive the tuple <int,int,double> from master?
      // gsInfo<<"Job on thread "<<my_rank<<" on level "<<l<<" with id "<<k<<" and value "<<val<<"\n";
      if (std::get<3>(tuple) < 0.75 && std::get<3>(tuple) > 0.25)
      {
        gsInfo<<"[MPI process "<<my_rank<<"] special!"<<"\n";
      }

      double maxSleep = 10;
      double sleepTime = maxSleep * next_random();

      sleep(sleepTime);
      gsInfo<<"slept for "<<sleepTime<<" seconds\n";

      gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to    "<<0<<": ID: "<<std::get<1>(tuple)<<"; x = "<<std::get<2>(tuple)<<"; y = "<<std::get<3>(tuple)<<"; level = "<<std::get<0>(tuple)<<" (slept for "<<sleepTime<<" seconds)\n";
      comm.send(&tuple,1,0,tag);

    }
  }


  return 0;
}