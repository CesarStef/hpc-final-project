/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */


#include "stencil_parallel.h"


// ------------------------------------------------------------------
// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  MPI_Comm myCOMM_WORLD;
  int  Rank, Ntasks;
  uint neighbours[4];

  int  Niterations;
  int  periodic;
  vec2_t S, N;
  
  int      Nsources;
  int      Nsources_local;
  vec2_t  *Sources_local;
  double   energy_per_source;

  plane_t   planes[2];  
  buffers_t buffers[2];
  
  int output_energy_stat_perstep;
  
  /* initialize MPI envrionment */
  {
    int level_obtained;
    
    // NOTE: change MPI_FUNNELED if appropriate
    //
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
    if ( level_obtained < MPI_THREAD_FUNNELED ) {
      // printf("MPI_thread level obtained is %d instead of %d\n",
	    //  level_obtained, MPI_THREAD_FUNNELED );
      MPI_Finalize();
      exit(1); 
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);
  }
  
  
  /* argument checking and setting */
  int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
			 neighbours, &Niterations,
			 &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
			 &planes[0], &buffers[0] );

  if ( ret )
  {
    // printf("task %d is opting out with termination code %d\n",
    //   Rank, ret );
    
    MPI_Finalize();
    return 0;
  }

  int current = OLD;
  double t1 = MPI_Wtime();   /* take wall-clock time */
  
  //printf("task %d starts the computation with %d sources, and will perform %d iterations\n", Rank, Nsources_local, Niterations);

  for (int iter = 0; iter < Niterations; ++iter)
  {
      MPI_Request reqs[8];
      
      /* new energy from sources */
      inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );


      /* -------------------------------------- */

      // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position
      
      // [A] Point buffer pointers to the correct positions in the plane data
      // The plane data layout is (size[_x_]+2) x (size[_y_]+2) with a halo frame
      // Row stride is (size[_x_]+2)
      uint sx = planes[current].size[_x_];
      uint sy = planes[current].size[_y_];
      uint stride = sx + 2;

      // NORTH halo: first inner row (row index 1), columns 1..sx
      buffers[SEND][NORTH] = &planes[current].data[1 * stride + 1];
      buffers[RECV][NORTH] = &planes[current].data[0 * stride + 1];

      // SOUTH halo: last inner row (row index sy), columns 1..sx
      buffers[SEND][SOUTH] = &planes[current].data[sy * stride + 1];
      buffers[RECV][SOUTH] = &planes[current].data[(sy + 1) * stride + 1];

      // EAST and WEST halos are not contiguous, so allocate/use separate buffers
      // Fill SEND buffer for WEST (column index 1)
      if ( neighbours[WEST] != MPI_PROC_NULL && buffers[SEND][WEST] != NULL )
        for ( uint j = 0; j < sy; j++ )
          buffers[SEND][WEST][j] = planes[current].data[(j + 1) * stride + 1];

      // Fill SEND buffer for EAST (column index sx)
      if ( neighbours[EAST] != MPI_PROC_NULL && buffers[SEND][EAST] != NULL )
        for ( uint j = 0; j < sy; j++ )
          buffers[SEND][EAST][j] = planes[current].data[(j + 1) * stride + sx];


      // [B] perfoem the halo communications
      // use Isend / Irecv
      //         --> can you overlap communication and compution in this way?

      int req_count = 0;

      // NORTH communications
      if ( neighbours[NORTH] != MPI_PROC_NULL ) {
        MPI_Isend( buffers[SEND][NORTH], sx, MPI_DOUBLE, neighbours[NORTH], NORTH, myCOMM_WORLD, &reqs[req_count++] );
        MPI_Irecv( buffers[RECV][NORTH], sx, MPI_DOUBLE, neighbours[NORTH], SOUTH, myCOMM_WORLD, &reqs[req_count++] );
      }

      // SOUTH communications
      if ( neighbours[SOUTH] != MPI_PROC_NULL ) {
        MPI_Isend( buffers[SEND][SOUTH], sx, MPI_DOUBLE, neighbours[SOUTH], SOUTH, myCOMM_WORLD, &reqs[req_count++] );
        MPI_Irecv( buffers[RECV][SOUTH], sx, MPI_DOUBLE, neighbours[SOUTH], NORTH, myCOMM_WORLD, &reqs[req_count++] );
      }

      // WEST communications
      if ( neighbours[WEST] != MPI_PROC_NULL ) {
        MPI_Isend( buffers[SEND][WEST], sy, MPI_DOUBLE, neighbours[WEST], WEST, myCOMM_WORLD, &reqs[req_count++] );
        MPI_Irecv( buffers[RECV][WEST], sy, MPI_DOUBLE, neighbours[WEST], EAST, myCOMM_WORLD, &reqs[req_count++] );
      }

      // EAST communications
      if ( neighbours[EAST] != MPI_PROC_NULL ) {
        MPI_Isend( buffers[SEND][EAST], sy, MPI_DOUBLE, neighbours[EAST], EAST, myCOMM_WORLD, &reqs[req_count++] );
        MPI_Irecv( buffers[RECV][EAST], sy, MPI_DOUBLE, neighbours[EAST], WEST, myCOMM_WORLD, &reqs[req_count++] );
      }

      //printf("task %d initiated %d communication requests for iteration %d\n", Rank, req_count, iter);

      MPI_Waitall( req_count, reqs, MPI_STATUSES_IGNORE );

      //printf("task %d completed all communications for iteration %d\n", Rank, iter);

      // [C] copy the halo data from RECV buffers back into the plane frame
      // WEST received data -> halo column 0
      if ( neighbours[WEST] != MPI_PROC_NULL )
        for ( uint j = 0; j < sy; j++ )
          planes[current].data[(j + 1) * stride + 0] = buffers[RECV][WEST][j];

      // EAST received data -> halo column sx+1
      if ( neighbours[EAST] != MPI_PROC_NULL )
        for ( uint j = 0; j < sy; j++ )
          planes[current].data[(j + 1) * stride + sx + 1] = buffers[RECV][EAST][j];
      
      // [C] copy the haloes data

      /* --------------------------------------  */
      /* update grid points */
      
      update_plane( periodic, N, &planes[current], &planes[!current] );
      
      /* swap plane indexes for the new iteration */
      current = !current;
      
    }
  
  t1 = MPI_Wtime() - t1;

  output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  
  if ( Rank == 0 ) {
    printf("Total Computational time: %f seconds\n", t1);
  }
  
  memory_release( buffers, planes );

  MPI_Finalize();

  // if ( Rank == 0 ) {
  //   printf("Total execution time: %f seconds\n", t1);
  // }

  return 0;
}


/* ==========================================================================
   =                                                                        =
   =   routines called within the integration loop                          =
   ========================================================================== */





/* ==========================================================================
   =                                                                        =
   =   initialization                                                       =
   ========================================================================== */


uint simple_factorization( uint, int *, uint ** );

int initialize_sources( int       ,
			int       ,
			MPI_Comm  *,
			uint      [2],
			int       ,
			int      *,
			vec2_t  ** );


int memory_allocate ( const int       *,
		      const vec2_t     N,
          buffers_t *,
          plane_t   *,
          int verbose);
		      

int initialize ( MPI_Comm *Comm,
		 int      Me,                  // the rank of the calling process
		 int      Ntasks,              // the total number of MPI ranks
		 int      argc,                // the argc from command line
		 char   **argv,                // the argv from command line
		 vec2_t  *S,                   // the size of the plane
		 vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
		 int     *periodic,            // periodic-boundary tag
		 int     *output_energy_stat,
		 int     *neighbours,          // four-int array that gives back the neighbours of the calling task
		 int     *Niterations,         // how many iterations
		 int     *Nsources,            // how many heat sources
		 int     *Nsources_local,
		 vec2_t **Sources_local,
		 double  *energy_per_source,   // how much heat per source
		 plane_t *planes,
		 buffers_t *buffers
		 )
{
  int halt = 0;
  int ret;
  int verbose = 0;
  
  // ··································································
  // set deffault values

  (*S)[_x_]         = 10000;
  (*S)[_y_]         = 10000;
  *periodic         = 0;
  *Nsources         = 4;
  *Nsources_local   = 0;
  *Sources_local    = NULL;
  *Niterations      = 1000;
  *energy_per_source = 1.0;

  if ( planes == NULL ) {
    // manage the situation
  }

  planes[OLD].size[0] = planes[OLD].size[0] = 0;
  planes[NEW].size[0] = planes[NEW].size[0] = 0;
  
  for ( int i = 0; i < 4; i++ )
    neighbours[i] = MPI_PROC_NULL;

  for ( int b = 0; b < 2; b++ )
    for ( int d = 0; d < 4; d++ )
      buffers[b][d] = NULL;
  
  // ··································································
  // process the commadn line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': (*S)[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': (*S)[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'o': *output_energy_stat = (atoi(optarg) > 0);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'v': verbose = atoi(optarg);
	    break;

	  case 'h': {
	    if ( Me == 0 )
	      printf( "\nvalid options are ( values btw [] are the default values ):\n"
		      "-x    x size of the plate [10000]\n"
		      "-y    y size of the plate [10000]\n"
		      "-e    how many energy sources on the plate [4]\n"
		      "-E    how many energy sources on the plate [1.0]\n"
		      "-n    how many iterations [1000]\n"
		      "-p    whether periodic boundaries applies  [0 = false]\n\n"
		      );
	    halt = 1; }
	    break;
	    
	    
	  case ':': printf( "option -%c requires an argument\n", optopt);
	    break;
	    
	  case '?': printf(" -------- help unavailable ----------\n");
	    break;
	  }
      }

    if ( opt == -1 )
      break;
  }

  if ( halt )
    return 1;
  
  
  // ··································································
  /*
   * here we should check for all the parms being meaningful
   *
   */

  // ...
  //TODO: add checks for the parameters
  
  // ··································································
  /*
   * find a suitable domain decomposition
   * very simple algorithm, you may want to
   * substitute it with a better one
   *
   * the plane Sx x Sy will be solved with a grid
   * of Nx x Ny MPI tasks
   */

  vec2_t Grid;
  double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
  int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) );

  
  if ( dimensions == 1 )
    {
      if ( (*S)[_x_] >= (*S)[_y_] )
	Grid[_x_] = Ntasks, Grid[_y_] = 1;
      else
	Grid[_x_] = 1, Grid[_y_] = Ntasks;
    }
  else
    {
      int   Nf;
      uint *factors;
      uint  first = 1;
      ret = simple_factorization( Ntasks, &Nf, &factors );
      
      for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ )
	first *= factors[i];

      if ( (*S)[_x_] > (*S)[_y_] )
	Grid[_x_] = Ntasks/first, Grid[_y_] = first;
      else
	Grid[_x_] = first, Grid[_y_] = Ntasks/first;
    }

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];
  

  // ··································································
  // my cooridnates in the grid of processors
  //
  int X = Me % Grid[_x_];
  int Y = Me / Grid[_x_];

  // ··································································
  // find my neighbours
  //

  if ( Grid[_x_] > 1 )
    {  
      if ( *periodic ) {       
	neighbours[EAST]  = Y*Grid[_x_] + (Me + 1 ) % Grid[_x_];
	neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
      
      else {
	neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
	neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
    }

  if ( Grid[_y_] > 1 )
  {
    if ( *periodic ) {      
      neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
      neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; 
    }
    else {    
      neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
      neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); 
    }
  }

  // ··································································
  // the size of my patch
  //

  /*
   * every MPI task determines the size sx x sy of its own domain
   * REMIND: the computational domain will be embedded into a frame
   *         that is (sx+2) x (sy+2)
   *         the outern frame will be used for halo communication or
   */
  
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_];
  uint r = (*S)[_x_] % Grid[_x_];
  mysize[_x_] = s + (X < r);
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r);

  planes[OLD].size[0] = mysize[0];
  planes[OLD].size[1] = mysize[1];
  planes[NEW].size[0] = mysize[0];
  planes[NEW].size[1] = mysize[1];
  

  if ( verbose > 0 )
    {
      if ( Me == 0 ) {
	      // printf("Tasks are decomposed in a grid %d x %d\n\n",
        // Grid[_x_], Grid[_y_] );
        fflush(stdout);
      }
      MPI_Barrier(*Comm);
      for ( int t = 0; t < Ntasks; t++ )
      {
        // if ( t == Me )
        //   {
        //     printf("Task %4d :: "
        //     "\tgrid coordinates : %3d, %3d\n"
        //     "\tneighbours: N %4d    E %4d    S %4d    W %4d\n",
        //     Me, X, Y,
        //     neighbours[NORTH], neighbours[EAST],
        //     neighbours[SOUTH], neighbours[WEST] );
        //     fflush(stdout);
        //   }

        MPI_Barrier(*Comm);
      }
    }

  
  // ··································································
  // allocae the needed memory
  //
  // if ( verbose > 0 ){
  //   printf("task %d is allocating memory for planes and buffers\n", Me);
  //   fflush(stdout);
  // }
  ret = memory_allocate( neighbours, mysize, buffers, planes, verbose );
  

  // ··································································
  // allocae the heat sources
  //
  // if ( verbose > 0 ){
  //   printf("task %d is initializing heat sources\n", Me);
  //   fflush(stdout);
  // }
  ret = initialize_sources( Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local );
  
  
  return 0;  
}


uint simple_factorization( uint A, int *Nfactors, uint **factors )
/*
 * rought factorization;
 * assumes that A is small, of the order of <~ 10^5 max,
 * since it represents the number of tasks
 #
 */
{
  int N = 0;
  int f = 2;
  uint _A_ = A;

  while ( f < A )
    {
      while( _A_ % f == 0 ) {
	N++;
	_A_ /= f; }

      f++;
    }

  *Nfactors = N;
  uint *_factors_ = (uint*)malloc( N * sizeof(uint) );

  N   = 0;
  f   = 2;
  _A_ = A;

  while ( f < A )
    {
      while( _A_ % f == 0 ) {
	_factors_[N++] = f;
	_A_ /= f; }
      f++;
    }

  *factors = _factors_;
  return 0;
}


int initialize_sources( int       Me,
			int       Ntasks,
			MPI_Comm *Comm,
			vec2_t    mysize,
			int       Nsources,
			int      *Nsources_local,
			vec2_t  **Sources )
{
  srand48(42 ^ Me); // different seed for each task to get different sources
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  
  if ( Me == 0 )
    {
      for ( int i = 0; i < Nsources; i++ )
	tasks_with_sources[i] = (int)lrand48() % Ntasks;
    }
  
  MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );

  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ )
    nlocal += (tasks_with_sources[i] == Me);
  *Nsources_local = nlocal;
  
  if ( nlocal > 0 )
    {
      vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
      for ( int s = 0; s < nlocal; s++ )
      {
        helper[s][_x_] = 1 + lrand48() % mysize[_x_];
        helper[s][_y_] = 1 + lrand48() % mysize[_y_];
      }
      *Sources = helper;
    }
  
  free( tasks_with_sources );

  // printf("task %d has %d sources\n", Me, nlocal);

  return 0;
}



int memory_allocate (
    const int *neighbours,
    const vec2_t     N,
    buffers_t *buffers_ptr,
    plane_t   *planes_ptr,
    int verbose
  )
{
    /*
      here you allocate the memory buffers that you need to
      (i)  hold the results of your computation
      (ii) communicate with your neighbours

      The memory layout that I propose to you is as follows:

      (i) --- calculations
      you need 2 memory regions: the "OLD" one that contains the
      results for the step (i-1)th, and the "NEW" one that will contain
      the updated results from the step ith.

      Then, the "NEW" will be treated as "OLD" and viceversa.

      These two memory regions are indexed by *plate_ptr:

      planew_ptr[0] ==> the "OLD" region
      plamew_ptr[1] ==> the "NEW" region


      (ii) --- communications

      you may need two buffers (one for sending and one for receiving)
      for each one of your neighnours, that are at most 4:
      north, south, east amd west.      

      To them you need to communicate at most mysizex or mysizey
      daouble data.

      These buffers are indexed by the buffer_ptr pointer so
      that

      (*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
      (*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
      
      --->> Of course you can change this layout as you prefer
      
     */

  if (planes_ptr == NULL ) {
    //printf("error: invalid pointer for planes\n");
    return 1;
  }

  if (buffers_ptr == NULL ) {
    //printf("error: invalid pointer for buffers\n");
    return 1;
  }

  if( verbose > 0 ) {
    //printf("task is allocating memory for planes and buffers\n");
    fflush(stdout);
  }
  // ··················································
  // allocate memory for data
  // we allocate the space needed for the plane plus a contour frame
  // that will contains data form neighbouring MPI tasks
  unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2);

  planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[OLD].data == NULL )
    {
      //printf("error allocating memory for plane data\n");
      return 1;
    }
  memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) );

  planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[NEW].data == NULL )
    {
      //printf("error allocating memory for plane data\n");
      return 1;
    }
  memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );

  // ··················································
  // buffers for north and south communication 
  // are not really needed
  //
  // in fact, they are already contiguous, just the
  // first and last line of every rank's plane
  //
  // you may just make some pointers pointing to the
  // correct positions
  //

  // or, if you preer, just go on and allocate buffers
  // also for north and south communications

  // ··················································
  // allocate buffers

  if ( neighbours[WEST] != MPI_PROC_NULL ) {
    buffers_ptr[SEND][WEST] = (double*)malloc( N[_y_] * sizeof(double) );
    buffers_ptr[RECV][WEST] = (double*)malloc( N[_y_] * sizeof(double) );
    memset( buffers_ptr[SEND][WEST], 0, N[_y_] * sizeof(double) );
    memset( buffers_ptr[RECV][WEST], 0, N[_y_] * sizeof(double) );
  }

  if ( neighbours[EAST] != MPI_PROC_NULL ) {
    buffers_ptr[SEND][EAST] = (double*)malloc( N[_y_] * sizeof(double) );
    buffers_ptr[RECV][EAST] = (double*)malloc( N[_y_] * sizeof(double) );
    memset( buffers_ptr[SEND][EAST], 0, N[_y_] * sizeof(double) );
    memset( buffers_ptr[RECV][EAST], 0, N[_y_] * sizeof(double) );
  }

  // if(verbose > 0){
  //   printf("task is done allocating memory for planes and buffers\n");
  //   fflush(stdout);
  // }

  return 0;
}

  int memory_release (buffers_t *buffers, plane_t   *planes)
  {
    if ( planes != NULL )
      {
        if ( planes[OLD].data != NULL )
          free (planes[OLD].data);
        if ( planes[NEW].data != NULL )
          free (planes[NEW].data);
      }
    return 0;
  }

  int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm )
{

  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  
  if ( Me == 0 )
    {
      printf( "total injected energy is %g, "
      "system energy is %g "
      "( in avg %g per grid point)\n",
      budget,
      tot_system_energy,
      tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
    }

  return 0;
}