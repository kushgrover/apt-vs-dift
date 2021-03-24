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
 * Sample uniformly over all the transitions
 */
public final class HeuristicNextStateUniform extends HeuristicNextState
{
	public HeuristicNextStateUniform(ModelGenerator pmg, StateUpdate stateUpdate)
	{
		super(pmg, stateUpdate);
	}

	@Override public State sample(State state, int action, int depth) throws PrismException
	{
		modelGenerator.exploreState(state);
		int trans = modelGenerator.getNumTransitions(action);
		int offset = random.nextInt(trans);

		updateChoice(state, action, offset);
		return modelGenerator.computeTransitionTarget(action, offset);
	}

	@Override public State sample(State state, int action) throws PrismException
	{
		return sample(state, action, 0);
	}

	@Override public NextState getType() {
		return NextState.UNIFORM;
	}
}
