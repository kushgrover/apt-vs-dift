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

package heuristics;

import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import parser.State;
import parser.Values;
import parser.VarList;
import parser.ast.RewardStruct;
import parser.type.Type;
import prism.ModelGenerator;
import prism.ModelType;
import prism.PrismException;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public final class CachedModelGenerator implements ModelGenerator
{
	private ModelGenerator modelGenerator = null;
	private final Map<State, CacheContainer> stateToCache = new HashMap<>();
	private State currentState = null;
	private CacheContainer currentCache = null;
	private final Map<State, State> mergedStateToRepresentativeMap = new HashMap<>();

	public CachedModelGenerator(ModelGenerator modelGenerator)
	{
		this.modelGenerator = modelGenerator;
	}

	public ModelGenerator getOriginalModelGenerator()
	{
		return modelGenerator;
	}

	@Override public boolean hasSingleInitialState()
	{
		throw new UnsupportedOperationException();
	}

	@Override public List<State> getInitialStates()
	{
		throw new UnsupportedOperationException();
	}

	@Override public void setSomeUndefinedConstants(Values someValues)
	{
		throw new UnsupportedOperationException();
	}

	@Override public State getInitialState() throws PrismException
	{
		return getStateRepresentative(modelGenerator.getInitialState());
	}

	@Override public Values getConstantValues()
	{
		throw new UnsupportedOperationException();
	}

	@Override public boolean containsUnboundedVariables()
	{
		throw new UnsupportedOperationException();
	}

	@Override public int getNumVars()
	{
		throw new UnsupportedOperationException();
	}

	@Override public void exploreState(State state) throws PrismException
	{
		assert state != null;
		final State newState = getStateRepresentative(state);
		if (newState.equals(currentState)) {
			assert currentCache == getFromCache(currentState);
			return;
		}
		currentState = newState;
		final CacheContainer storedCache = getFromCache(currentState);
		if (storedCache == null) {
			currentCache = buildCache();
			putInCache(currentState, currentCache);
		} else {
			currentCache = storedCache;
		}
	}

	@Override public List<String> getVarNames()
	{
		throw new UnsupportedOperationException();
	}

	@Override public List<Type> getVarTypes()
	{
		throw new UnsupportedOperationException();
	}

	@Override public State getExploreState()
	{
		return currentState;
	}

	@Override public int getVarIndex(String name)
	{
		throw new UnsupportedOperationException();
	}

	@Override public String getVarName(int i)
	{
		throw new UnsupportedOperationException();
	}

	@Override public int getNumLabels()
	{
		throw new UnsupportedOperationException();
	}

	@Override public List<String> getLabelNames()
	{
		throw new UnsupportedOperationException();
	}

	@Override public String getLabelName(int i)
	{
		throw new UnsupportedOperationException();
	}

	@Override public int getLabelIndex(String label)
	{
		throw new UnsupportedOperationException();
	}

	@Override public int getNumRewardStructs()
	{
		throw new UnsupportedOperationException();
	}

	@Override public List<String> getRewardStructNames()
	{
		throw new UnsupportedOperationException();
	}

	@Override public int getRewardStructIndex(String name)
	{
		throw new UnsupportedOperationException();
	}

	@Override public RewardStruct getRewardStruct(int i)
	{
		throw new UnsupportedOperationException();
	}

	@Override public boolean rewardStructHasTransitionRewards(int i)
	{
		throw new UnsupportedOperationException();
	}

	@Override public VarList createVarList()
	{
		throw new UnsupportedOperationException();
	}

	@Override public Object getChoiceAction(int i)
	{
		if (currentCache == null) {
			throw new IllegalStateException();
		}
		return currentCache.choiceToAction.get(i);
	}

	@Override public boolean isLabelTrue(String label)
	{
		throw new UnsupportedOperationException();
	}

	@Override public boolean isLabelTrue(int i)
	{
		throw new UnsupportedOperationException();
	}

	@Override public double getStateReward(int r, State state)
	{
		throw new UnsupportedOperationException();
	}

	@Override public double getStateActionReward(int r, State state, Object action)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public int getNumChoices()
	{
		assert currentCache != null && currentCache == getFromCache(currentState);
		return currentCache.numChoices;
	}

	@Override
	public int getNumTransitions(int i)
	{
		assert currentCache != null && currentCache == getFromCache(currentState);
		return currentCache.numTrans.getInt(i);
	}

	@Override
	public double getTransitionProbability(int i, int offset)
	{
		assert currentCache != null && currentCache == getFromCache(currentState);
		return currentCache.choiceToDistribution.get(i).get(offset);
	}

	@Override
	public State computeTransitionTarget(int i, int offset)
	{
		assert currentCache != null && currentCache == getFromCache(currentState);
		return currentCache.choiceToSuccessor.get(i).get(offset);
	}

	@Override
	public int getNumTransitions()
	{
		assert currentCache != null && currentCache == getFromCache(currentState);
		return currentCache.numTransitions;
	}

	public int getPlayerForState()
	{
		assert currentCache != null && currentCache == getFromCache(currentState);
		return currentCache.player;
	}

	@Override
	public ModelType getModelType()
	{
		throw new UnsupportedOperationException();
		/* See commit 6f5a1531c3a96528311af060e47c11d18121cf6f for what was here */
	}

	@Override
	public String getTransitionAction(int i, int offset)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public String getTransitionAction(int i)
	{
		throw new UnsupportedOperationException();
	}

	public void updateCache(State state, List<Object2DoubleMap<State>> actions)
	{
		CacheContainer cache = getFromCache(state);
		if (cache == null) {
			throw new IllegalStateException("Cannot update empty cache");
		}
		putInCache(state, buildCache(actions));
	}

	public State mergeStates(Set<State> states, List<Object2DoubleMap<State>> actions) throws PrismException
	{
		if (states.isEmpty()) {
			throw new PrismException("Can't merge empty set");
		}
		final Iterator<State> stateIterator = states.iterator();
		final State representative = getStateRepresentative(stateIterator.next());
		assert representative == getStateRepresentative(representative) : representative + " is part of a MEC and not the representative";
		assert stateToCache.containsKey(representative) : representative + " not in cache";

		final CacheContainer representativeCache = buildCache(actions);
		putInCache(representative, representativeCache);

		// TODO Add a cached reverse lookup
		Set<State> transitiveStates = new HashSet<>(states);
		for (Map.Entry<State, State> entry : mergedStateToRepresentativeMap.entrySet()) {
			if (transitiveStates.contains(entry.getValue())) {
				transitiveStates.add(entry.getKey());
			}
		}

		if (transitiveStates.size() > 1) {
			mergedStateToRepresentativeMap.put(representative, representative);
			for (State state : transitiveStates) {
				if (state == representative) {
					continue;
				}
				stateToCache.remove(state);
				if (state.equals(currentState)) {
					currentState = representative;
					currentCache = representativeCache;
				}
				mergedStateToRepresentativeMap.put(state, representative);
			}
		}

		return representative;
	}

	private CacheContainer buildCache() throws PrismException
	{
		CacheContainer cache = new CacheContainer();

		modelGenerator.exploreState(currentState);
		cache.numChoices = modelGenerator.getNumChoices();
		cache.numTrans = new IntArrayList();
		cache.choiceToDistribution = new Int2ObjectArrayMap<>();
		cache.choiceToSuccessor = new Int2ObjectArrayMap<>();
		cache.choiceToAction = new Int2ObjectArrayMap<>();

		if (cache.numChoices > 0) {
			for (int i = 0; i < cache.numChoices; i++) {
				cache.numTrans.add(modelGenerator.getNumTransitions(i));

				Int2DoubleMap dist = new Int2DoubleArrayMap();
				cache.choiceToDistribution.put(i, dist);

				Int2ObjectArrayMap<State> successor = new Int2ObjectArrayMap<>();
				cache.choiceToSuccessor.put(i, successor);

				cache.choiceToAction.put(i, modelGenerator.getChoiceAction(i));

				for (int j = 0; j < cache.numTrans.getInt(i); j++) {
					dist.put(j, modelGenerator.getTransitionProbability(i, j));
					successor.put(j, modelGenerator.computeTransitionTarget(i, j));
				}
			}
			cache.numTransitions = modelGenerator.getNumTransitions();
		} else {
			//Deadlock
			/*cache.numTrans = new int[1];
			cache.numTrans[0] = 1;
			Map<Integer,Double> dist = new HashMap<Integer,Double>();
			cache.choiceToDistribution.put(0, dist);
			Map<Integer, State> succ = new HashMap<Integer, State>();
			cache.choiceToSuccessor.put(0, succ);
			dist.put(0, 1.0);
			succ.put(0, currentState);
			cache.numTransitions = 1;*/
			// TODO: Do we really need the deadlock handling?
			cache.numTransitions = 0;
		}
		if (modelGenerator.getModelType() == ModelType.SMG) {
			// TODO: PORT modelExplorer.getPlayerForState().
			// cache.player = getPlayerForState();
		}
		return cache;
	}

	private CacheContainer buildCache(List<Object2DoubleMap<State>> actions)
	{
		CacheContainer cache = new CacheContainer();

		cache.numChoices = actions.size();
		cache.numTrans = new IntArrayList();
		cache.choiceToDistribution = new Int2ObjectArrayMap<>();
		cache.choiceToSuccessor = new Int2ObjectArrayMap<>();
		cache.choiceToAction = new Int2ObjectArrayMap<>();

		int trans = 0;
		for (int i = 0; i < cache.numChoices; i++) {
			cache.numTrans.add(actions.get(i).size());

			Int2DoubleMap dist = new Int2DoubleArrayMap();
			cache.choiceToDistribution.put(i, dist);

			Int2ObjectArrayMap<State> succ = new Int2ObjectArrayMap<>();
			cache.choiceToSuccessor.put(i, succ);

			Iterator<Object2DoubleMap.Entry<State>> it = actions.get(i).object2DoubleEntrySet().iterator();
			int j = 0;

			while (it.hasNext()) {
				Object2DoubleMap.Entry<State> e = it.next();
				final State state = getStateRepresentative(e.getKey());
				dist.put(j, e.getDoubleValue());
				succ.put(j, state);
				j++;
			}
			trans += cache.numTrans.getInt(i);
		}
		cache.numTransitions = trans;
		if (modelGenerator.getModelType() == ModelType.SMG) {
			throw new IllegalStateException("Cannot update cache for SMGs");
			//cache.player = modelExplorer.getPlayerForState();
		}
		return cache;
	}

	private CacheContainer getFromCache(State state)
	{
		return stateToCache.get(getStateRepresentative(state));
	}

	public State getStateRepresentative(State state) {
		State representative = mergedStateToRepresentativeMap.get(state);
		if (representative == null) {
			return state;
		}
		return representative;
	}

	private void putInCache(State state, CacheContainer cacheContainer)
	{
		assert state == getStateRepresentative(state);
		if (state.equals(currentState)) {
			currentCache = cacheContainer;
		}
		stateToCache.put(state, cacheContainer);
	}

	static final class CacheContainer
	{
		int numTransitions = -1;
		int numChoices = -1;
		IntList numTrans = null;
		int player = -1;

		Int2ObjectMap<Int2DoubleMap> choiceToDistribution = null;
		Int2ObjectMap<Int2ObjectMap<State>> choiceToSuccessor = null;
		Int2ObjectMap<Object> choiceToAction = null;
	}

}
