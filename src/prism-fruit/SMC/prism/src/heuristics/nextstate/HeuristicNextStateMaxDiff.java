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

package heuristics.nextstate;

import heuristics.nextstate.HeuristicNextStateFactory.NextState;
import heuristics.search.StateValue;
import heuristics.update.StateUpdate;
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

public final class HeuristicNextStateMaxDiff extends HeuristicNextState
{
	public HeuristicNextStateMaxDiff(ModelGenerator modelGenerator, StateUpdate stateUpdate)
	{
		super(modelGenerator, stateUpdate);
	}

	@Override public State sample(State state, int action, int depth) throws PrismException
	{
		modelGenerator.exploreState(state);
		int trans = modelGenerator.getNumTransitions(action);
		if (trans > 1) {
			double[] diffS = new double[trans];
			for (int i = 0; i < trans; i++) {
				State t = modelGenerator.computeTransitionTarget(action, i);
				Double p = modelGenerator.getTransitionProbability(action, i);
				StateValue sv;
				if (depth != -1) {
					sv = stateUpdate.getQValue(t, depth + 1);
				} else {
					sv = stateUpdate.getQValue(t);
				}
				double lowerBound = 0d;
				double upperBound = 1d;
				if (sv != null) {
					lowerBound = sv.lowerBound;
					upperBound = sv.upperBound;
				}
				diffS[i] = p * (upperBound - lowerBound);
			}
			return sampleFromDist(state, action, diffS);
		} else {
			updateChoice(state, action, 0);
			return modelGenerator.computeTransitionTarget(action, 0);
		}
	}

	@Override public State sample(State state, int action) throws PrismException
	{
		return sample(state, action, -1);
	}

	@Override public NextState getType()
	{
		return NextState.MAX_DIFF;
	}

}
