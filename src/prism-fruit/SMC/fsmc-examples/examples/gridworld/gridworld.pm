dtmc

const int  N; // grid size
const double p;

formula lfree = (janitor_col > 1) &(robot_col != janitor_col-1 | robot_row !=janitor_row);
formula rfree = (janitor_col < N) & (robot_col != janitor_col+1 | robot_row !=janitor_row);
formula ufree = (janitor_row > 1) &  (robot_row != janitor_row-1 | robot_col !=janitor_col);
formula dfree = (janitor_row < N) & (robot_row != janitor_row+1 | robot_col !=janitor_col);

module robot

robot_row   :[1..N] init N;
robot_col   :[1..N] init 1;

[move]  connection & (robot_col < N) & (janitor_col != robot_col+1 | robot_row !=janitor_row) -> 
	(robot_col'=robot_col+1); // first move right

[move]  connection & (robot_row > 1) & (robot_col = N) & (janitor_row != robot_row-1 | robot_col !=janitor_col) -> 
	(robot_row'=robot_row-1); // second move up

[move] connection & (robot_row = 1) & (robot_col = N) -> true; 

[move] !connection -> true;
endmodule


module janitor

janitor_row   :[1..N] init 1;
janitor_col   :[1..N] init N;


[move] lfree & rfree & ufree & dfree ->
	0.25 : (janitor_col'=janitor_col-1) +
	0.25 : (janitor_col'=janitor_col+1) +
	0.25 : (janitor_row'=janitor_row-1) +
	0.25 : (janitor_row'=janitor_row+1);

[move] !lfree & rfree & ufree & dfree ->
	1/3 : (janitor_col'=janitor_col+1) +
	1/3 : (janitor_row'=janitor_row-1) +
	1/3 : (janitor_row'=janitor_row+1);

[move] lfree & !rfree & ufree & dfree ->
	1/3 : (janitor_col'=janitor_col-1) +
	1/3 : (janitor_row'=janitor_row-1) +
	1/3 : (janitor_row'=janitor_row+1);

[move] lfree & rfree & !ufree & dfree ->
	1/3 : (janitor_col'=janitor_col-1) +
	1/3 : (janitor_col'=janitor_col+1) +
	1/3 : (janitor_row'=janitor_row+1);

[move] lfree & rfree & ufree & !dfree ->
	1/3 : (janitor_col'=janitor_col-1) +
	1/3 : (janitor_col'=janitor_col+1) +
	1/3 : (janitor_row'=janitor_row-1);

[move] !lfree & !rfree & ufree & dfree ->
	0.5 : (janitor_row'=janitor_row-1) +
	0.5 : (janitor_row'=janitor_row+1);

[move]  !lfree & rfree & ufree & !dfree ->
	0.5 : (janitor_col'=janitor_col+1) +
	0.5 : (janitor_row'=janitor_row-1);

[move] lfree & rfree & !ufree & !dfree  ->
	0.5 : (janitor_col'=janitor_col-1) +
	0.5 : (janitor_col'=janitor_col+1);

[move] lfree & !rfree & !ufree & dfree ->
	0.5 : (janitor_col'=janitor_col-1) +
	0.5 : (janitor_row'=janitor_row+1);

[move] lfree & !rfree & ufree & !dfree ->
	0.5 : (janitor_col'=janitor_col-1) +
	0.5 : (janitor_row'=janitor_row-1);

[move] !lfree & rfree & !ufree & dfree ->
	0.5 : (janitor_col'=janitor_col+1) +
	0.5 : (janitor_row'=janitor_row+1);

[move] !lfree & !rfree & !ufree & dfree ->
	(janitor_row'=janitor_row+1);

[move] !lfree & !rfree & ufree & !dfree ->
	(janitor_row'=janitor_row-1);

[move] !lfree & rfree & !ufree & !dfree ->
	(janitor_col'=janitor_col+1);

[move] lfree & !rfree & !ufree & !dfree ->
	(janitor_col'=janitor_col-1);

endmodule


module station 
connection : bool init true; // is the robot connected to the communication station?

	[move] connection -> p : (connection'=true) + (1-p) : (connection'=false);
        [move] !connection -> (connection'=false);

endmodule


rewards 
	(robot_row = 1) & (robot_col = N) : 0.9;
endrewards