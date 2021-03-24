//==============================================================================
//
//	Copyright (c) 2016-
//	Authors:
//	* Joachim Klein <klein@tcs.inf.tu-dresden.de> (TU Dresden)
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
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSet;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;
import strat.MDStrategy;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntConsumer;

/**
 * This is a minimal implementation of an explicitly stored labeled transition system.
 * Each target state of an edge is modeled as a choice, with a single transition for
 * this choice.
 */
public class LTSExplicit extends ModelExplicit implements LTS
{
	protected List<IntList> successors = new ArrayList<>();
	protected int numTransitions = 0;
	protected int maxNumChoices = 0;

	public LTSExplicit()
	{
		initialise(0);
	}

	public int addState()
	{
		successors.add(new IntArrayList());
		return numStates++;
	}

	public void addEdge(int s, int t)
	{
		successors.get(s).add(t);
		numTransitions++;
		maxNumChoices = Math.max(getNumChoices(s), maxNumChoices);
	}

	@Override
	public int getNumChoices(int s)
	{
		// one choice per successor for s
		return successors.get(s).size();
	}

	@Override
	public int getMaxNumChoices()
	{
		return maxNumChoices;
	}

	@Override
	public int getNumChoices()
	{
		return numTransitions;
	}

	@Override
	public Object getAction(int s, int i)
	{
		// we don't have action labels
		return null;
	}

	@Override
	public boolean areAllChoiceActionsUnique()
	{
		// as we don't assign action labels, they are not unique
		return false;
	}

	@Override
	public int getNumTransitions(int s, int i)
	{
		if (i < getNumChoices(s)) {
			// one transition per choice
			return 1;
		}
		throw new IllegalArgumentException();
	}

	public int getNumTransitions(int s)
	{
		return getNumChoices(s);
	}

	@Override
	public boolean allSuccessorsInSet(int s, int i, IntSet set)
	{
		// single successor for s, i
		return set.contains(successors.get(s).getInt(i));
	}

	@Override
	public boolean someSuccessorsInSet(int s, int i, IntSet set)
	{
		// single successor for s, i
		return set.contains(successors.get(s).getInt(i));
	}

	@Override
	public IntIterator getSuccessorsIterator(int s, int i)
	{
		// single successor for s, i
		return IntIterators.singleton(successors.get(s).getInt(i));
	}

	@Override
	public Model constructInducedModel(MDStrategy strat)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFileWithStrat(PrismLog out, IntSet mark, MDStrategy strat)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public IntIterator getSuccessorsIterator(int s)
	{
		return IntIterators.unmodifiable(successors.get(s).iterator());
	}

	@Override
	public void findDeadlocks(boolean fix) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public ModelType getModelType()
	{
		return ModelType.LTS;
	}

	@Override
	public int getNumTransitions()
	{
		return numTransitions;
	}

	@Override
	public void exportToPrismExplicitTra(PrismLog out)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	protected void exportTransitionsToDotFile(int s, PrismLog out)
	{
		getSuccessorsIterator(s).forEachRemaining((IntConsumer) successor -> out.println(s + " -> " + successor + ";"));
	}

	@Override
	public void exportToPrismLanguage(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}
}
