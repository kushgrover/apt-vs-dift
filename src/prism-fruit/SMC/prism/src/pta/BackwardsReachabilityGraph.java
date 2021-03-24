//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

package pta;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.Distribution;
import explicit.MDP;
import explicit.MDPSimple;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import prism.PrismException;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntConsumer;

public class BackwardsReachabilityGraph
{
	public List<LocZone> states;
	private IntList initialStates;
	private NatBitSet target;
	private List<List<List<IntList>>> trans;

	public class Edge
	{
		int index; // Edge index in transition
		int dest; // Destination location

		public Edge(int index, int dest)
		{
			this.index = index;
			this.dest = dest;
		}

		@Override
		public String toString()
		{
			return index + "/" + dest;
		}

		@Override
		public boolean equals(Object o)
		{
			return o instanceof Edge && ((Edge) o).index == index && ((Edge) o).dest == dest;
		}
	}

	public BackwardsReachabilityGraph()
	{
		initialStates = new IntArrayList();
		target = NatBitSets.set();
		trans = new ArrayList<>();
	}

	public void addState(List<Transition> trs)
	{
		int numTransitions = trs.size();
		List<List<IntList>> list = new ArrayList<>(numTransitions);
		for (int i = 0; i < numTransitions; i++) {
			int numEdges = trs.get(i).getNumEdges();
			List<IntList> list2 = new ArrayList<>(numEdges);
			for (int j = 0; j < numEdges; j++) {
				list2.add(new IntArrayList());
			}
			list.add(list2);
		}
		trans.add(list);
	}

	public void addInitialState(int s)
	{
		initialStates.add(s);
	}

	public void addTargetState(int s)
	{
		target.set(s);
	}

	public IntList getInitialStates()
	{
		return initialStates;
	}

	public NatBitSet getTarget()
	{
		return target;
	}

	public void addTransition(int src, int tr, int i, int dest)
	{
		IntList list = trans.get(src).get(tr).get(i);
		if (!list.contains(dest))
			list.add(dest);
	}

	public List<List<IntList>> getList(int src)
	{
		return trans.get(src);
	}

	public MDP buildMDP(PTA pta)
	{
		MDPSimple mdp = new MDPSimple(states.size() + 1); // + sink
		int src = -1;
		for (List<List<IntList>> list : trans) {
			src++;
			int tr = -1;
			for (List<IntList> list2 : list) {
				tr++;
				Distribution distr = new Distribution();
				double prob, rest = 0;
				int j = -1;
				for (IntList dests : list2) {
					j++;
					prob = pta.getTransitions(states.get(src).loc).get(tr).getEdges().get(j).getProbability();
					if (dests.size() > 1) {
						int sNew = mdp.addState();
						distr.add(sNew, prob);
						dests.forEach((IntConsumer) dest -> mdp.addChoice(sNew, new Distribution(dest, 1.0)));
					} else if (dests.size() == 1) {
						distr.add(dests.getInt(0), prob);
					} else {
						rest += prob;
					}
				}
				//if (rest > 0)
				//	distr.add(mdp.getNumStates() - 1, rest);
				mdp.addChoice(src, distr);
			}
		}
		// Add initial states
		initialStates.forEach((IntConsumer) mdp::addInitialState);
		// fix sink
		try {
			mdp.findDeadlocks(true);
		} catch (PrismException e) {
			// Never happens for MDPSimple
		}
		//log.println(mdp);
		return mdp;
	}

	public MDP buildMdpExpo(PTA pta)
	{
		MDPSimple mdp = new MDPSimple(states.size() + 1); // + sink
		int src = -1;
		for (List<List<IntList>> list : trans) {
			src++;
			int tr = -1;
			for (List<IntList> list2 : list) {
				tr++;
				int dests[] = new int[list2.size()];
				int size = 1;
				for (IntList aList2 : list2) {
					if (!aList2.isEmpty()) {
						size *= aList2.size();
					}
				}
				if (size > 6) {
					System.out.println(size + "!");
					System.out.println(list2);
					for (IntList list3 : list2) {
						list3.forEach((IntConsumer) x -> System.out.println(x + ":" + states.get(x)));
					}

				}
				buildMdpExpo(mdp, pta, src, tr, list2, 0, dests);
			}
		}
		// Add initial states
		initialStates.forEach((IntConsumer) mdp::addInitialState);
		// fix sink
		try {
			mdp.findDeadlocks(true);
		} catch (PrismException e) {
			// Never happens for MDPSimple
		}
		//log.println(mdp);
		return mdp;
	}

	public void buildMdpExpo(MDPSimple mdp, PTA pta, int src, int tr, List<IntList> list2, int i, int dests[])
	{
		if (i == dests.length) {
			//log.print(src + "/" + tr);
			Distribution distr = new Distribution();
			double prob, rest = 0;
			if (dests.length > 0) {
				for (int j = 0; j < dests.length; j++) {
					prob = pta.getTransitions(states.get(src).loc).get(tr).getEdges().get(j).getProbability();
					if (!list2.get(j).isEmpty()) {
						//log.print(" " + prob + ":" + dests[j]);
						distr.add(dests[j], prob);
					} else {
						rest += prob;
					}
				}
			}
			//if (rest > 0)
			//	distr.add(mdp.getNumStates() - 1, rest);
			mdp.addChoice(src, distr);
			//log.println();
		} else {
			IntList list3 = list2.get(i);
			if (list3.isEmpty()) {
				buildMdpExpo(mdp, pta, src, tr, list2, i + 1, dests);
			} else {
				list3.forEach((IntConsumer) dest -> {
					dests[i] = dest;
					buildMdpExpo(mdp, pta, src, tr, list2, i + 1, dests);
				});
			}
		}
	}

	@Override
	public String toString()
	{
		return trans.toString();
	}
}
