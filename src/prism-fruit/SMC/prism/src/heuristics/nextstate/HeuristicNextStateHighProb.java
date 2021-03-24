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
import heuristics.update.StateUpdate;
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

/**
 * Sample according to the distribution
 */
public final class HeuristicNextStateHighProb extends HeuristicNextState
{
	public HeuristicNextStateHighProb(ModelGenerator pmg, StateUpdate stateUpdate)
	{
		super(pmg, stateUpdate);
	}

	@Override public State sample(State state, int action, int depth) throws PrismException
	{
		// TODO: Handle deadlocks
		modelGenerator.exploreState(state);
		int trans = modelGenerator.getNumTransitions(action);
		double[] probabilities = new double[trans];
		for (int i = 0; i < trans; i++) {
			probabilities[i] = modelGenerator.getTransitionProbability(action, i);
		}
		return sampleFromDist(state, action, probabilities);
	}

	@Override public State sample(State state, int action) throws PrismException
	{
		return sample(state, action, -1);
	}

	@Override public NextState getType()
	{
		return NextState.HIGH_PROB;
	}

}
