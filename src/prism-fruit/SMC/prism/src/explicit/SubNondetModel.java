//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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

package explicit;

import common.FastUtils;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.Int2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import parser.State;
import parser.Values;
import parser.VarList;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;
import strat.MDStrategy;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.IntConsumer;

/**
 * Class for creating a sub-model of any NondetModel, please note the translate* methods
 * used to translate between state ids for model and sub-model. Created sub-model will have new
 * state numbering from 0 to number of states in the sub model.
 */
public class SubNondetModel implements NondetModel
{
	/**
	 * (Optionally) the stored predecessor relation. Becomes inaccurate after the model is changed!
	 */
	protected PredecessorRelation predecessorRelation;
	private NondetModel model = null;
	private NatBitSet states = null;
	private Int2ObjectMap<NatBitSet> actions = null;
	private NatBitSet initialStates = null;
	private List<State> statesList = null;
	private Int2IntMap stateLookupTable = new Int2IntOpenHashMap();
	private Int2ObjectMap<Int2IntMap> actionLookupTable = new Int2ObjectOpenHashMap<>();
	private Int2IntMap inverseStateLookupTable = new Int2IntOpenHashMap();
	private int numTransitions = 0;
	private int maxNumChoices = 0;
	private int numChoices = 0;

	public SubNondetModel(NondetModel model, NatBitSet states, Int2ObjectMap<NatBitSet> actions, NatBitSet initialStates)
	{
		this.model = model;
		this.states = states;
		this.actions = actions;
		this.initialStates = initialStates;

		generateStatistics();
		generateLookupTable(states, actions);
	}

	@Override
	public ModelType getModelType()
	{
		return model.getModelType();
	}

	@Override
	public int getNumStates()
	{
		return states.size();
	}

	@Override
	public int getNumInitialStates()
	{
		return initialStates.size();
	}

	@Override
	public IntIterable getInitialStates()
	{
		IntList is = new IntArrayList();
		initialStates.forEach((IntConsumer) i -> is.add(translateState(i)));
		return is;
	}

	@Override
	public int getFirstInitialState()
	{
		return translateState(initialStates.firstInt());
	}

	@Override
	public boolean isInitialState(int i)
	{
		return initialStates.contains(translateState(i));
	}

	@Override
	public int getNumDeadlockStates()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public IntIterable getDeadlockStates()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public StateValues getDeadlockStatesList()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public List<State> getStatesList()
	{
		// We use lazy generation because in many cases the state list is not needed
		if (statesList == null) {
			statesList = generateSubStateList(states);
		}
		return statesList;
	}

	private List<State> generateSubStateList(NatBitSet states)
	{
		List<State> statesList = new ArrayList<>();
		IntIterator iterator = FastUtils.iterator(states, model.getNumStates());
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			statesList.add(model.getStatesList().get(i));
		}
		return statesList;
	}

	@Override
	public VarList getVarList()
	{
		// we can return the varList of the model, as we do not change
		// the variables in the model
		return model.getVarList();
	}

	@Override
	public Values getConstantValues()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public Set<String> getLabels()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public NatBitSet getLabelStates(String name)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean hasLabel(String name)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public int getNumTransitions()
	{
		return numTransitions;
	}

	@Override
	public IntIterator getSuccessorsIterator(int s)
	{
		int translateState = translateState(s);
		IntSet succs = new IntOpenHashSet();

		IntIterator iterator = FastUtils.iterator(actions.get(translateState), model.getNumChoices(translateState));
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			model.getSuccessorsIterator(translateState, i).forEachRemaining((int j) -> succs.add(inverseTranslateState(j)));
		}

		return succs.iterator();
	}

	@Override
	public boolean isSuccessor(int s1, int s2)
	{
		s1 = translateState(s1);
		s2 = translateState(s2);
		return model.isSuccessor(s1, s2);
	}

	@Override public void findDeadlocks(boolean fix) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicit(String baseFilename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicitTra(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicitTra(File file) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicitTra(PrismLog log)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFile(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFile(String filename, IntSet mark) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFile(PrismLog out)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFile(PrismLog out, IntSet mark)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFile(PrismLog out, IntSet mark, boolean showStates)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFileWithStrat(PrismLog out, IntSet mark, MDStrategy strat)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismLanguage(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportStates(int exportType, VarList varList, PrismLog log) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public String infoString()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public String infoStringTable()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public int getNumChoices(int s)
	{
		s = translateState(s);
		return actions.get(s).size();
	}

	@Override
	public int getMaxNumChoices()
	{
		return maxNumChoices;
	}

	@Override
	public int getNumChoices()
	{
		return numChoices;
	}

	@Override
	public Object getAction(int s, int i)
	{
		int sOriginal = translateState(s);
		int iOriginal = translateAction(s, i);

		return model.getAction(sOriginal, iOriginal);
	}

	@Override
	public boolean areAllChoiceActionsUnique()
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public int getNumTransitions(int s, int i)
	{
		int sOriginal = translateState(s);
		int iOriginal = translateAction(s, i);
		return model.getNumTransitions(sOriginal, iOriginal);
	}

	@Override
	public IntIterator getSuccessorsIterator(int s, int i)
	{
		int sOriginal = translateState(s);
		int iOriginal = translateAction(s, i);
		IntSet set = new IntOpenHashSet();
		model.getSuccessorsIterator(sOriginal, iOriginal).forEachRemaining((int j) -> set.add(inverseTranslateState(j)));
		return set.iterator();
	}

	private NatBitSet translateSet(NatBitSet set)
	{
		NatBitSet translatedNatBitSet = NatBitSets.setWithExpectedSize(set.size());
		set.forEach((IntConsumer) i -> translatedNatBitSet.set(translateState(i)));
		return translatedNatBitSet;
	}

	private void generateStatistics()
	{
		IntIterator iterator = FastUtils.iterator(states, model.getNumStates());
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			numTransitions += getTransitions(i);
			numChoices += actions.get(i).size();
			maxNumChoices = Math.max(maxNumChoices, model.getNumChoices(i));
		}
	}

	private int getTransitions(int state)
	{
		int transitions = 0;
		IntIterator iterator = FastUtils.iterator(actions.get(state), model.getNumChoices(state));
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			transitions += model.getNumTransitions(state, i);
		}
		return transitions;
	}

	@Override
	public Model constructInducedModel(MDStrategy strat)
	{
		throw new RuntimeException("Not implemented");
	}

	private void generateLookupTable(NatBitSet states, Int2ObjectMap<NatBitSet> actions)
	{
		IntIterator iterator = FastUtils.iterator(states, model.getNumStates());
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			inverseStateLookupTable.put(i, stateLookupTable.size());
			stateLookupTable.put(stateLookupTable.size(), i);
			Int2IntMap r = new Int2IntLinkedOpenHashMap();
			IntIterator availableActions = FastUtils.iterator(actions.get(i), model.getNumChoices(i));
			availableActions.forEachRemaining((IntConsumer) j -> r.put(r.size(), j));
			actionLookupTable.put(actionLookupTable.size(), r);
		}
	}

	public int translateState(int s)
	{
		return stateLookupTable.get(s);
	}

	private int inverseTranslateState(int s)
	{
		return inverseStateLookupTable.get(s);
	}

	public int translateAction(int s, int i)
	{
		return actionLookupTable.get(s).get(i);
	}

	@Override
	public boolean hasStoredPredecessorRelation()
	{
		return (predecessorRelation != null);
	}

	@Override
	public PredecessorRelation getPredecessorRelation(prism.PrismComponent parent, boolean storeIfNew)
	{
		if (predecessorRelation != null) {
			return predecessorRelation;
		}

		PredecessorRelation pre = PredecessorRelations.forModel(parent, this);

		if (storeIfNew) {
			predecessorRelation = pre;
		}
		return pre;
	}

	@Override
	public void clearPredecessorRelation()
	{
		predecessorRelation = null;
	}
}
