//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package heuristics.search;

import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.update.StateUpdate;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Stack;

import explicit.ModelExplorer;

import parser.State;
import prism.PrismException;

//Labelled version of LRTDP

public class HeuristicLRTDP extends Heuristic{

	private Map<State, BitSet> state2Solved = new HashMap<State, BitSet>();
	private Stack<VisitedState> visited = new Stack<VisitedState>();
	
	public HeuristicLRTDP(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
		super(mc, su, nextState, pme, min);
	}
	
	@Override
	public void heuristicStep(State s) throws PrismException {
		int depth = 0;
		visited.clear();
		while(!isSolved(s,depth)) {
			visited.push(new VisitedState(s, depth));
			int bestAction = stateUpdate.update(s, depth);
			s = nextState.sample(s,bestAction, depth);
			depth++;
			trialSteps++;
		}
		stateUpdate.setTargetOrBoundValue(s, depth);
		while(!visited.empty()) {
			VisitedState v = visited.pop();
			if(!checkSolved(v.state,v.depth)) {
				stateUpdate.update(v.state, v.depth);
				break;
			} else {
				stateUpdate.update(v.state, v.depth);
			}
		}
		if(trials % 10000 == 0) {
			reportProgress(trials, trialSteps);
		}
		trials++;
	}

	@Override
	public boolean isDone() throws PrismException {
		return checkSolved(initialState, 0);
	}
	
	private boolean checkSolved(State s, int depth) throws PrismException{
		if(isSolved(s, depth)) {
			markSolved(s, depth);
			return true;
		}
		double lowerBound = stateUpdate.getQValue(s,depth).getLowerBound();
		double upperBound = stateUpdate.getQValue(s,depth).getUpperBound();
		if(stateUpdate.same(lowerBound, upperBound)) {
			markSolved(s, depth);
			return true;
		}
		boolean solved = true;
		pme.queryState(s);
		int choices = pme.getNumChoices();
		for(int i=0;i<choices;i++) {
			pme.queryState(s);
			int numTrans = pme.getNumTransitions(i);
			for(int j=0;j<numTrans;j++) {
				pme.queryState(s);
				State succ = pme.computeTransitionTarget(i, j);
				solved = solved & isSolved(succ, depth+1);
			}
		}
		
		if(solved) {
			markSolved(s, depth);
			return true;
		} else {
			return false;
		}
	}
	
	private void markSolved(State s, int depth) throws PrismException{
		BitSet solvedAtDepth = state2Solved.get(s);
		if(solvedAtDepth == null) {
			solvedAtDepth = new BitSet();
			state2Solved.put(s, solvedAtDepth);
		}
		solvedAtDepth.set(depth, true);
	}
	
	private boolean isSolved(State s, int depth) throws PrismException{
		if(stateUpdate.isTarget(s) || (depth >= bound) || stateUpdate.isSelfLoop(s)) {
			markSolved(s, depth);
			return true;
		}
		if(state2Solved.get(s) == null) {
			return false;
		} else {
			return state2Solved.get(s).get(depth);
		}
	}
	
	class VisitedState {
		State state;
		int depth;
		
		public VisitedState(State s, int d) {
			this.state = s;
			this.depth = d;
		}
		
		public String toString() {
			return "State " + state + " depth " + depth;
 		}
	}
}
