dtmc

module m2
  s : [0..4] init 0;

  [] s<4 -> 1/3 : (s'=0) + 1/3 : (s'=s) + 1/3 : (s'=s+1);
  [] s=4 -> 1 : (s'=s);
endmodule
