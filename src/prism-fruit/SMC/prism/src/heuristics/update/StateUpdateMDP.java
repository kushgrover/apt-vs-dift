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

package heuristics.update;

import heuristics.CachedModelGenerator;
import heuristics.search.StateValue;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import parser.State;
import prism.PrismException;

import java.util.Iterator;
import java.util.Random;

public final class StateUpdateMDP extends StateUpdate
{
	private static final Random random = new Random();

	public StateUpdateMDP(CachedModelGenerator modelGenerator, StateValueContainer container, int bound, boolean min, double epsilon) throws PrismException
	{
		super(modelGenerator, container, bound, min, epsilon);
	}

	/**
	 * Updates the Q value for state state and returns an
	 * (epsilon-close) choice which will give this value
	 * @param state
	 * @param depth
	 * @return
	 * @throws PrismException
	 */
	@Override public int update(State state, int depth) throws PrismException
	{
		modelGenerator.exploreState(state);
		int numChoices = modelGenerator.getNumChoices();
		if (numChoices == 0) {
			return -1;
		}

		int bestAction = 0;
		StateValue bestBounds = getActionValueBounds(state, 0, depth + 1);
		for (int i = 1; i < numChoices; i++) {
			StateValue successorBounds = getActionValueBounds(state, i, depth + 1);
			if (min) {
				if (successorBounds.lowerBound < bestBounds.lowerBound) {
					bestAction = i;
					bestBounds = successorBounds;
				}
			} else {
				// Pick the action that gives the largest upper bound
				if (successorBounds.upperBound > bestBounds.upperBound) {
					bestAction = i;
					bestBounds = successorBounds;
				}
			}
		}

		setQValue(state, bestBounds, depth);
		return bestAction;
	}

	/**
	 * Updates the Q value for state state and returns an
	 * (epsilon-close) choice which will give this value
	 * @param state
	 * @return
	 * @throws PrismException
	 */
	@Override public int update(State state) throws PrismException
	{
		modelGenerator.exploreState(state);
		int numChoices = modelGenerator.getNumChoices();
		if (numChoices == 0) {
			return -1;
		}

		// Look for the best action with respect to the corresponding successor bounds. If min == true, search for the action which gives the lowest lower
		// bound, similarly if min == false, find the action that gives the largest upper bound.
		// Furthermore, keep track of all actions which are similar (i.e. bounds are only epsilon far apart)

		Int2ObjectMap<StateValue> similarBounds = new Int2ObjectArrayMap<>();
		int bestChoice = 0;
		StateValue bestBounds = getActionValueBounds(state, bestChoice);
		for (int choice = 1; choice < numChoices; choice++) {
			StateValue successorBounds = getActionValueBounds(state, choice);
			if (min) {
				if (successorBounds.lowerBound < bestBounds.lowerBound) {
					if (same(successorBounds.lowerBound, bestBounds.lowerBound)) {
						similarBounds.put(bestChoice, bestBounds);
					} else {
						similarBounds.clear();
					}
					bestChoice = choice;
					bestBounds = successorBounds;
				}
			} else {
				if (successorBounds.upperBound > bestBounds.upperBound) {
					if (same(successorBounds.upperBound, bestBounds.upperBound)) {
						similarBounds.put(bestChoice, bestBounds);
					} else {
						similarBounds.clear();
					}
					bestChoice = choice;
					bestBounds = successorBounds;
				}
			}
		}
		assert !similarBounds.containsKey(bestChoice);

		// similarBounds only gathered all which are potentially close - filter out the invalid ones.
		final Iterator<Int2ObjectMap.Entry<StateValue>> iterator = similarBounds.int2ObjectEntrySet().iterator();
		while (iterator.hasNext()) {
			final Int2ObjectMap.Entry<StateValue> entry = iterator.next();
			final StateValue actionBounds = entry.getValue();
			if (min) {
				assert actionBounds.lowerBound >= bestBounds.lowerBound;
				if (!same(actionBounds.lowerBound, bestBounds.lowerBound)) {
					iterator.remove();
				}
			} else {
				assert actionBounds.upperBound <= bestBounds.upperBound;
				if (!same(actionBounds.upperBound, bestBounds.upperBound)) {
					iterator.remove();
				}
			}
		}

		// Now select which action to use if there are similar ones
		assert !similarBounds.containsKey(bestChoice);
		if (!similarBounds.isEmpty()) {
			int randomActionIndex = random.nextInt(similarBounds.size() + 1);
			if (randomActionIndex != 0) {
				// 0 means "take the current bestAction", 1 is "take the first from the similar bounds set" etc.
				final Iterator<Int2ObjectMap.Entry<StateValue>> boundsIterator = similarBounds.int2ObjectEntrySet().iterator();
				while (randomActionIndex > 1) {
					boundsIterator.next();
					randomActionIndex--;
				}
				final Int2ObjectMap.Entry<StateValue> entry = boundsIterator.next();
				bestChoice = entry.getIntKey();
				bestBounds = entry.getValue();
			}
		}

		setQValue(state, bestBounds);
		return bestChoice;
	}
}
