//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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

package explicit;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntSet;
import prism.PrismLog;
import strat.MDStrategy;

import java.util.function.IntConsumer;
import java.util.stream.IntStream;

/**
 * Interface for (abstract) classes that provide (read-only) access to an explicit-state model with nondeterminism.
 */
public interface NondetModel extends Model
{
	// Accessors

	default void forEachTransition(int s, int i, IntConsumer consumer)
	{
		getSuccessors(s, i).forEach(consumer);
	}

	/**
	 * Get the number of nondeterministic choices in state s.
	 */
	int getNumChoices(int s);

	/**
	 * Get the maximum number of nondeterministic choices in any state.
	 */
	int getMaxNumChoices();

	/**
	 * Get the total number of nondeterministic choices over all states.
	 */
	int getNumChoices();

	/**
	 * Get the action label (if any) for choice {@code i} of state {@code s}.
	 */
	Object getAction(int s, int i);

	/**
	 * Do all choices in in each state have a unique action label?
	 */
	boolean areAllChoiceActionsUnique();

	/**
	 * Get the number of transitions from choice {@code i} of state {@code s}.
	 */
	int getNumTransitions(int s, int i);

	/**
	 * Check if all the successor states from choice {@code i} of state {@code s} are in the set {@code set}.
	 * @param s The state to check
	 * @param i Choice index
	 * @param set The set to test for inclusion
	 */
	default boolean allSuccessorsInSet(int s, int i, IntSet set) {
		return IntIterators.all(getSuccessorsIterator(s, i), set::contains);
	}

	/**
	 * Check if some successor state from choice {@code i} of state {@code s} is in the set {@code set}.
	 * @param s The state to check
	 * @param i Choice index
	 * @param set The set to test for inclusion
	 */
	default boolean someSuccessorsInSet(int s, int i, IntSet set)  {
		return IntIterators.any(getSuccessorsIterator(s, i), set::contains);
	}

	/**
	 * Get an iterator over the successor states from choice {@code i} of state {@code s}.
	 * @param s The state
	 * @param i Choice index
	 */
	IntIterator getSuccessorsIterator(int s, int i);

	default IntStream getSuccessorsIntStream(int s, int i) {
		IntStream.Builder successorBuilder = IntStream.builder();
		getSuccessorsIterator(s, i).forEachRemaining((IntConsumer) successorBuilder::add);
		return successorBuilder.build();
	}

	default IntCollection getSuccessors(int s, int i) {
		return new IntArrayList(getSuccessorsIterator(s, i));
	}

	/**
	 * Construct a model that is induced by applying strategy {@code strat} to this model.
	 * Note that the "new" model may be just an implicit (read-only) representation.
	 * @param strat (Memoryless) strategy to use
	 */
	Model constructInducedModel(MDStrategy strat);

	/**
	 * Export to a dot file, highlighting states in 'mark' and choices for a (memoryless) strategy.
	 */
	void exportToDotFileWithStrat(PrismLog out, IntSet mark, MDStrategy strat);
}