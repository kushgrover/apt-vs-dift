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

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.update.StateUpdate;
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.Map;

//Labelled version of LRTDP

public class HeuristicLrtdp extends Heuristic
{
	private long modelCheckingTime = 0;
	private final Map<State, NatBitSet> stateToSolved = new HashMap<>();
	private final Deque<VisitedState> visited = new ArrayDeque<>();

	public HeuristicLrtdp(HeuristicsSMGModelChecker smgModelChecker, StateUpdate stateUpdate, HeuristicNextState nextState, ModelGenerator modelGenerator,
			boolean min) throws PrismException
	{
		super(smgModelChecker, stateUpdate, nextState, modelGenerator, min);
	}

	@Override
	protected void heuristicStart() throws PrismException
	{
		modelCheckingTime = System.currentTimeMillis();
	}

	@Override
	protected void heuristicStop() throws PrismException
	{
		long duration = System.currentTimeMillis() - modelCheckingTime;
		smgModelChecker.getLog().println();
		smgModelChecker.getLog().println("Heuristic model checking time in " + ((double)duration/1000)  + " secs.");
	}

	@Override
	public void exploreAndUpdate(State s) throws PrismException
	{
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
		while (!visited.isEmpty()) {
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
	public boolean isDone() throws PrismException
	{
		return checkSolved(initialState, 0);
	}

	private boolean checkSolved(State s, int depth) throws PrismException
	{
		if(isSolved(s, depth)) {
			markSolved(s, depth);
			return true;
		}
		double lowerBound = stateUpdate.getQValue(s, depth).lowerBound;
		double upperBound = stateUpdate.getQValue(s, depth).upperBound;
		if(stateUpdate.same(lowerBound, upperBound)) {
			markSolved(s, depth);
			return true;
		}
		boolean solved = true;
		modelGenerator.exploreState(s);
		int choices = modelGenerator.getNumChoices();
		for(int i=0;i<choices;i++) {
			modelGenerator.exploreState(s);
			int numTrans = modelGenerator.getNumTransitions(i);
			for(int j=0;j<numTrans;j++) {
				modelGenerator.exploreState(s);
				State succ = modelGenerator.computeTransitionTarget(i, j);
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

	private void markSolved(State s, int depth) throws PrismException
	{
		NatBitSet solvedAtDepth = stateToSolved.get(s);
		if(solvedAtDepth == null) {
			solvedAtDepth = NatBitSets.set();
			stateToSolved.put(s, solvedAtDepth);
		}
		solvedAtDepth.set(depth, true);
	}

	private boolean isSolved(State s, int depth) throws PrismException
	{
		if(stateUpdate.isTarget(s) || (depth >= bound) || stateUpdate.isSelfLoop(s)) {
			markSolved(s, depth);
			return true;
		}
		if(stateToSolved.get(s) == null) {
			return false;
		} else {
			return stateToSolved.get(s).contains(depth);
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
