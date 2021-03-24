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

import heuristics.update.StateUpdate;
import prism.ModelGenerator;
import prism.PrismLangException;

public final class HeuristicNextStateFactory
{
	private HeuristicNextStateFactory() {}

	public enum NextState {
		UNIFORM, HIGH_PROB, DETERMINISTIC, MAX_DIFF
	}

	public static HeuristicNextState getHeuristicNextState(NextState nextState, ModelGenerator modelGenerator, StateUpdate stateUpdate)
			throws PrismLangException
	{
		switch(nextState) {
			case UNIFORM:
				return new HeuristicNextStateUniform(modelGenerator, stateUpdate);
			case HIGH_PROB:
				return new HeuristicNextStateHighProb(modelGenerator, stateUpdate);
			case DETERMINISTIC:
				return new HeuristicNextStateDeterministic(modelGenerator, stateUpdate);
			case MAX_DIFF:
				return new HeuristicNextStateMaxDiff(modelGenerator, stateUpdate);
		}
		throw new PrismLangException("Unknown next state heuristic");
	}
}
