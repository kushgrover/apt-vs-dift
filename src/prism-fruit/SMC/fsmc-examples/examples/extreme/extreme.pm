// synchronous leader election protocol  (itai & Rodeh)
// dxp/gxn 25/01/01

dtmc

// CONSTANTS
const int N; // number of processes
const double p=0.5;

module extreme

	s : [0..N];

	[] s<=N-1 -> (1-p) : (s'=0) + p : (s'=s+1);
	[] s=N   -> (s'=0);

endmodule

label "goal"= (s=N);

rewards
	true : 1;
endrewards
