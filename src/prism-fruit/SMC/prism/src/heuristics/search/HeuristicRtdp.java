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
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

import java.util.ArrayDeque;
import java.util.Deque;

public class HeuristicRtdp extends Heuristic
{
	private long modelCheckingTime = 0;
	private final Deque<State> visited = new ArrayDeque<>();


	public HeuristicRtdp(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelGenerator pmg, boolean min) throws PrismException
	{
		super(mc, su, nextState, pmg, min);
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
	public void exploreAndUpdate(State state) throws PrismException
	{
		int depth = 0;
		while(!stateUpdate.isTargetState(state) && depth < bound) {
			visited.push(state);
			int bestAction = stateUpdate.update(state, depth);
			state = nextState.sample(state, bestAction, depth);
			depth++;
			trialSteps++;
		}
		stateUpdate.setTargetOrBoundValue(state, depth);
		while (!visited.isEmpty()) {
			State visitedState = visited.pop();
			assert depth > 0;
			depth -= 1;
			stateUpdate.update(visitedState, depth);
		}
		if(trials % 10000 == 0) {
			reportProgress(trials, trialSteps);
		}
		trials++;
	}

	@Override
	public boolean isDone() throws PrismException
	{
		int depth = getBound() == -1 ? -1 : 0;
		StateValue sv = stateUpdate.getQValue(initialState, depth);
		if(sv != null) {
			double lowerBoundInitialState = sv.lowerBound;
			double upperBoundInitialState = sv.upperBound;
			return stateUpdate.same(lowerBoundInitialState, upperBoundInitialState);
		}
		return false;
	}

}