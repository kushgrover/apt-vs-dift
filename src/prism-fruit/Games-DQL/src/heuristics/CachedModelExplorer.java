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

import explicit.ModelExplorer;
import parser.State;
import prism.ModelType;
import prism.PrismException;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class CachedModelExplorer implements ModelExplorer{

	private ModelExplorer modelExplorer = null;
	private Map<State, CacheContainer> state2Cache = new HashMap<>();
	private State currentState = null;

	public CachedModelExplorer(ModelExplorer me) {
		modelExplorer = me;
	}

	@Override
	public State getDefaultInitialState() throws PrismException
	{
		return modelExplorer.getDefaultInitialState();
	}

	@Override
	public void queryState(State state) throws PrismException
	{
		currentState = state;
		CacheContainer cache = getFromCache(state);
		if(cache == null) {
			modelExplorer.queryState(state);
			putInCache(state, buildCache());
		}
	}

	@Override
	public int getNumChoices() throws PrismException
	{
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.numChoices;
		} else {
			throw new IllegalStateException();
		}
	}

	@Override
	public int getNumTransitions(int i) throws PrismException
	{
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.numTrans[i];
		} else {
			throw new IllegalStateException();
		}

	}

	@Override
	public double getTransitionProbability(int i, int offset) throws PrismException
	{
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.choice2Distribution.get(i).get(offset);
		} else {
			throw new IllegalStateException();
		}

	}

	@Override
	public State computeTransitionTarget(int i, int offset) throws PrismException
	{
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.choice2Successor.get(i).get(offset);
		} else {
			throw new IllegalStateException();
		}

	}

	@Override
	public int getNumTransitions() throws PrismException
	{
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.numTransitions;
		} else {
			throw new IllegalStateException();
		}
	}

	public int getPlayerForState() throws PrismException {
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.player;
		} else {
			throw new IllegalStateException();
		}
	}

	public Map<String, Integer> getPlayerMapping() throws PrismException{
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.players;
		} else {
			throw new IllegalStateException();
		}
	}

	public ModelType getModelType() throws PrismException {
		CacheContainer cache = getFromCache(currentState);
		if(cache != null) {
			return cache.modelType;
		} else {
			throw new IllegalStateException();
		}
	}

	@Override
	public void queryState(State state, double time) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public double getTransitionProbability(int i) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	//Body added by Maxi on 15.01.18; if there is only one successor, return it; else throw exception
	@Override
	public State computeTransitionTarget(int i) throws PrismException
	{
		CacheContainer cache = getFromCache(currentState);
		if(cache != null && cache.choice2Successor.get(i).size()==1) {
			return cache.choice2Successor.get(i).get(0);
		} else {
			throw new PrismException("More than one successor, but no offset specified");
		}
	}

	@Override
	public String getTransitionAction(int i, int offset) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public String getTransitionAction(int i) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

//	@Override
//	public ModelExplorer clone() {
//		return new CachedModelExplorer(modelExplorer.clone());
//	}

	public void updateCache(State s, List<Map<State, Double>> actions) throws PrismException{
		CacheContainer c = state2Cache.get(s);
		if(c != null) {
			putInCache(s, buildCache(c, actions));
		} else {
			throw new IllegalStateException("Cannot update empty cache");
		}
	}

	private CacheContainer buildCache() throws PrismException{
		CacheContainer cache = new CacheContainer();
		cache.numChoices = modelExplorer.getNumChoices();
		cache.numTrans = new int[modelExplorer.getNumChoices()];
		cache.choice2Distribution = new HashMap<Integer, Map<Integer,Double>>();
		cache.choice2Successor = new HashMap<Integer, Map<Integer,State>>();

		if(cache.numChoices > 0) {
			for(int i=0;i<cache.numChoices;i++) {
				cache.numTrans[i] = modelExplorer.getNumTransitions(i);
				Map<Integer,Double> dist = new HashMap<Integer,Double>();
				cache.choice2Distribution.put(i, dist);
				Map<Integer,State> succ = new HashMap<Integer,State>();
				cache.choice2Successor.put(i, succ);

				for(int j=0;j<cache.numTrans[i];j++) {
					dist.put(j, modelExplorer.getTransitionProbability(i, j));
					succ.put(j, modelExplorer.computeTransitionTarget(i, j));
				}
			}
			cache.numTransitions = modelExplorer.getNumTransitions();
		} else {
			//Deadlock
			cache.numTrans = new int[1];
			cache.numTrans[0] = 1;
			Map<Integer,Double> dist = new HashMap<Integer,Double>();
			cache.choice2Distribution.put(0, dist);
			Map<Integer,State> succ = new HashMap<Integer,State>();
			cache.choice2Successor.put(0, succ);
			dist.put(0, 1.0);
			succ.put(0, currentState);
			cache.numTransitions = 1;
		}
		if(modelExplorer.getModelType() == ModelType.SMG) {
			cache.player = modelExplorer.getPlayerForState();
			cache.players = modelExplorer.getPlayerMapping();
		}
		if(modelExplorer.getModelType() == ModelType.MDP) {
			cache.player = 1;
			Map<String, Integer> playerMapping = new HashMap<String, Integer>();
			playerMapping.put(String.valueOf(1), 1); //TODO: 12.01.18 Put 1 here, earlier was SMG.PLAYER_1; no idea if this will break
			cache.players = playerMapping;
		}
		return cache;
	}

	private CacheContainer buildCache(List<Map<State, Double>> actions) throws PrismException{
		CacheContainer cache = new CacheContainer();
		cache.numChoices = actions.size();
		cache.numTrans = new int[actions.size()];
		cache.choice2Distribution = new HashMap<Integer, Map<Integer,Double>>();
		cache.choice2Successor = new HashMap<Integer, Map<Integer,State>>();
		int trans = 0;
		for(int i=0;i<cache.numChoices;i++) {
			cache.numTrans[i] = actions.get(i).size();

			Map<Integer,Double> dist = new HashMap<Integer,Double>();
			cache.choice2Distribution.put(i, dist);
			Map<Integer,State> succ = new HashMap<Integer,State>();
			cache.choice2Successor.put(i, succ);

			Iterator<Entry<State, Double>> it = actions.get(i).entrySet().iterator();
			int j = 0;
			while(it.hasNext()) {
				Entry<State, Double> e = it.next();
				dist.put(j, e.getValue());
				succ.put(j, e.getKey());
				j++;
			}
			trans += cache.numTrans[i];
		}
		cache.numTransitions = trans;
		if(modelExplorer.getModelType() == ModelType.SMG) {
			throw new IllegalStateException("Cannot update cache for SMGs");
			//cache.player = modelExplorer.getPlayerForState();
		}
		return cache;
	}

	private CacheContainer buildCache(CacheContainer oldCache, List<Map<State, Double>> actions) throws PrismException{
		CacheContainer cache = new CacheContainer();
		cache.numChoices = actions.size();
		cache.numTrans = new int[actions.size()];
		cache.choice2Distribution = new HashMap<Integer, Map<Integer,Double>>();
		cache.choice2Successor = new HashMap<Integer, Map<Integer,State>>();
		int trans = 0;
		cache.modelType = oldCache.modelType;

		for(int i=0;i<cache.numChoices;i++) {
			cache.numTrans[i] = actions.get(i).size();

			Map<Integer,Double> dist = new HashMap<Integer,Double>();
			cache.choice2Distribution.put(i, dist);
			Map<Integer,State> succ = new HashMap<Integer,State>();
			cache.choice2Successor.put(i, succ);

			Iterator<Entry<State, Double>> it = actions.get(i).entrySet().iterator();
			int j = 0;
			while(it.hasNext()) {
				Entry<State, Double> e = it.next();
				dist.put(j, e.getValue());
				succ.put(j, e.getKey());
				j++;
			}
			trans += cache.numTrans[i];
		}
		cache.numTransitions = trans;

		cache.player = oldCache.player;
		cache.players = oldCache.players;

		return cache;
	}


	private CacheContainer getFromCache(State s) {
		return state2Cache.get(s);
	}

	private void putInCache(State s, CacheContainer c) {
		state2Cache.put(s, c);
	}

	class CacheContainer {
		int numTransitions = -1;
		int numChoices = -1;
		int numTrans[] = null;
		int player = -1;
		ModelType modelType;
		Map<String, Integer> players;

		Map<Integer, Map<Integer,Double>> choice2Distribution = null;
		Map<Integer, Map<Integer,State>> choice2Successor = null;
	}

	public State getSucc(int a){
		return (State) state2Cache.get(currentState).choice2Successor.get(a).get(0);
	}

}
