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
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

import java.util.HashMap;
import java.util.Map;

/**
 * Round robin sampling
 */
public final class HeuristicNextStateDeterministic extends HeuristicNextState
{
	private final Map<State, Int2ObjectMap<Int2IntMap>> stateToLastTransBounded = new HashMap<>();

	public HeuristicNextStateDeterministic(ModelGenerator modelGenerator, StateUpdate stateUpdate)
	{
		super(modelGenerator, stateUpdate);
	}

	@Override public State sample(State state, int action, int depth) throws PrismException
	{
		modelGenerator.exploreState(state);
		Int2IntMap map = stateToLastTransBounded.computeIfAbsent(state, s -> new Int2ObjectOpenHashMap<>())
				.computeIfAbsent(depth, d -> new Int2IntOpenHashMap());
		int transitionOffset = 0;
		if (map.containsKey(action)) {
			transitionOffset = map.get(action);
			transitionOffset++;
			int numTrans = modelGenerator.getNumTransitions(action);
			if (transitionOffset >= numTrans) {
				transitionOffset = 0;
			}
		}
		map.put(action, transitionOffset);
		updateChoice(state, action, transitionOffset);
		return modelGenerator.computeTransitionTarget(action, transitionOffset);
	}

	/**
	 * Chooses transitions from action one after the other
	 *
	 * @param state
	 * @param action
	 * @return
	 * @throws PrismException
	 */
	@Override public State sample(State state, int action) throws PrismException
	{
		modelGenerator.exploreState(state);
		int transitionOffset = 0;
		Int2IntMap lastTransitionMap = stateToLastChosenTrans.get(state);
		if (lastTransitionMap != null) {
			if (lastTransitionMap.containsKey(action)) {
				int lastTransition = lastTransitionMap.get(action);
				transitionOffset = (lastTransition + 1) % modelGenerator.getNumTransitions(action);
			}
		} else {
			lastTransitionMap = new Int2IntOpenHashMap();
			stateToLastChosenTrans.put(state, lastTransitionMap);
		}
		lastTransitionMap.put(action, transitionOffset);

		updateChoice(state, action, transitionOffset);
		return modelGenerator.computeTransitionTarget(action, transitionOffset);
	}

	@Override public NextState getType()
	{
		return NextState.DETERMINISTIC;
	}
}
