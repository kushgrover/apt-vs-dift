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

import java.util.LinkedList;
import java.util.Stack;

import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.search.HeuristicLRTDP.VisitedState;
import heuristics.update.StateUpdate;

import explicit.ModelExplorer;

import parser.State;
import prism.PrismException;

public class HeuristicILAO extends Heuristic{

	private LinkedList<State> fringe = new LinkedList<State>();
	
	public HeuristicILAO(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
		super(mc, su, nextState, pme, min);
	}
	
	public void compute() throws PrismException{
		State initial = pme.getDefaultInitialState();
		double lowerBoundInitialState = 0;
		double upperBoundInitialState = 1;
		boolean done = false;
		
		stateUpdate.setQValue(initial, new StateValue(lowerBoundInitialState, upperBoundInitialState));
		fringe.add(initial);
		
		heuristicStart();
		
		while(!done) {
			State s = fringe.removeLast();
			Stack<State> visited = dfs(s);
			
			while(!visited.empty()) {
				State v = visited.pop();
				//TODO: If I ever look at this again, need to adjust this to stateUpdate changes stateUpdate.update(v);
			}
			
			done = isDone();
		}
		reportProgress(trials, trialSteps);
		heuristicStop();
	}
	
	private Stack<State> dfs(State s) {
		return null;
	}
	
	@Override
	public void heuristicStep(State s) throws PrismException {
		int depth = 0;
		while(!stateUpdate.isTarget(s) && !stateUpdate.isSelfLoop(s) && depth < bound) {
			int bestAction = stateUpdate.update(s, depth);
			s = nextState.sample(s,bestAction, depth);
			trialSteps++;
		}
		stateUpdate.setTargetOrBoundValue(s, depth);
		if(trials % 10000 == 0) {
			reportProgress(trials, trialSteps);
		}
		trials++;
	}

	@Override
	public boolean isDone() throws PrismException {
		StateValue sv = stateUpdate.getQValue(initialState);
		if(sv != null) {
			double lowerBoundInitialState = sv.getLowerBound();
			double upperBoundInitialState = sv.getUpperBound();
			return stateUpdate.same(lowerBoundInitialState, upperBoundInitialState);
		}
		return false;
	}

}