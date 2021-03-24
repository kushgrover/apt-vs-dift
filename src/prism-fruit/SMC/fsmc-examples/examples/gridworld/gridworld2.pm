dtmc

const int  N; // grid size
const double p;



module robot

robot_row   :[1..N] init N;
robot_col   :[1..N] init 1;

[move]  (robot_col < N) & (janitor_col != robot_col+1 | robot_row !=janitor_row) -> 
	(robot_col'=robot_col+1); // first move right

[move]  (robot_row > 1) & (robot_col = N) & (janitor_row != robot_row-1 | robot_col !=janitor_col) -> 
	(robot_row'=robot_row-1); // second move up

[move] (robot_row = 1) & (robot_col = N) -> true; 
endmodule


module janitor

janitor_row   :[1..N] init 1;
janitor_col   :[1..N] init N;


[move]  (janitor_row > 1) & (robot_row != janitor_row-1 | robot_col !=janitor_col) -> 
	(janitor_row'=janitor_row-1); // move up

[move]  (janitor_row < N) & (robot_row != janitor_row+1 | robot_col !=janitor_col) -> 
	(janitor_row'=janitor_row+1); // move down

[move]  (janitor_col > 1) & (robot_col != janitor_col-1 | robot_row !=janitor_row) -> 
	(janitor_col'=janitor_col-1); // move left


[move]  (janitor_col < N) & (robot_col != janitor_col+1 | robot_row !=janitor_row) -> 
	(janitor_col'=janitor_col+1); // move right
endmodule


module station 
connection : bool init true; // is the robot connected to the communication station?

	[move] connection -> p : (connection'=true) + (1-p) : (connection'=false);
        [move] !connection -> (connection'=false);

endmodule