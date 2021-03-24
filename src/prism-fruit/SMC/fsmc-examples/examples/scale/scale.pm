dtmc

// CONSTANTS
const int N; // number of processes
    const int M;
const double pmin=0.5;

module extreme

	s : [0..N+M]; // tail 0..N, BSCC N+1..N+M
	b : [0..1];

	[] s> 0 & s <= N -> (1-pmin) : (s'=s) + pmin : (s'=s+1);
	[] s = N+1 -> 	0.5: (s'=s+1) + 0.5 : (s'=N+M);
	[] s = N+M ->   0.5: (s'=N+1) + 0.5 : (s'=N+M-1);
	[] s>N+1 & s <N+M -> 0.5 : (s'=s+1) + 0.5 : (s'=s-1);

	[] s= 0 -> 1/2 : (s'=1) & (b'=0) + 1/2 : (s'=1) & (b'=1);
endmodule

label "goal"= s>N & b=1;

rewards
	s>N: 1;
endrewards
