/*
 * Copyright (C) 2016 (Salomon Sickert)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package owl;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrays;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import owl.automaton.Automaton;
import owl.automaton.acceptance.GeneralizedBuchiAcceptance;
import owl.automaton.acceptance.NoneAcceptance;
import owl.automaton.edge.Edge;
import owl.automaton.ldba.LimitDeterministicAutomaton;
import owl.ltl.EquivalenceClass;
import owl.ltl.LabelledFormula;
import owl.run.env.DefaultEnvironment;
import owl.translations.Optimisation;
import owl.translations.ltl2ldba.LTL2LDBAFunction;
import owl.translations.ltl2ldba.breakpointfree.FGObligations;
import owl.translations.ltl2ldba.breakpointfree.GeneralizedBreakpointFreeState;
import parser.ast.Expression;
import prism.PrismLog;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

public class LDBAWrapper
{
	private static final Function<LabelledFormula, LimitDeterministicAutomaton<EquivalenceClass, GeneralizedBreakpointFreeState, GeneralizedBuchiAcceptance,
			FGObligations>> TRANSLATOR;

	static {
		EnumSet<Optimisation> optimisations = EnumSet.allOf(Optimisation.class);
		optimisations.remove(Optimisation.REMOVE_EPSILON_TRANSITIONS);
		TRANSLATOR = LTL2LDBAFunction.createGeneralizedBreakpointFreeLDBABuilder(new DefaultEnvironment(false, false), optimisations);
	}

	private final LimitDeterministicAutomaton<EquivalenceClass, GeneralizedBreakpointFreeState, GeneralizedBuchiAcceptance, FGObligations> ldba;
	private final Object2IntMap<Object> numbering;
	private final List<Object> map;
	private final Int2ObjectMap<int[]> jumpCache;
	private final int smallestAcceptingState;

	private final Automaton<EquivalenceClass, NoneAcceptance> initialComponent;
	private final Automaton<GeneralizedBreakpointFreeState, GeneralizedBuchiAcceptance> acceptingComponent;

	public LDBAWrapper(LimitDeterministicAutomaton<EquivalenceClass, GeneralizedBreakpointFreeState, GeneralizedBuchiAcceptance, FGObligations> ldba)
	{
		this.ldba = ldba;

		// The non-accepting trap state is stored as -1.
		numbering = new Object2IntOpenHashMap<>(ldba.size());
		numbering.defaultReturnValue(-2);
		map = new ArrayList<>(ldba.size());

		initialComponent = ldba.getInitialComponent();
		acceptingComponent = ldba.getAcceptingComponent();

		// Register initial state as state 0;
		if (initialComponent.getInitialStates().isEmpty()) {
			get(acceptingComponent.getInitialState());
		} else {
			get(initialComponent.getInitialState());
		}

		ldba.getInitialComponent().getStates().forEach(this::get);
		smallestAcceptingState = map.size();
		ldba.getAcceptingComponent().getStates().forEach(this::get);
		jumpCache = new Int2ObjectOpenHashMap<>();
	}

	public static LDBAWrapper create(Expression ltl, PrismLog logger)
	{
		logger.println("\nBuilding deterministic automaton (for " + ltl + ")...");

		LabelledFormula ltlFormula = ExpressionToFormula.convert(ltl);
		logger.println("Got owl-formula " + ltlFormula);
		long time = System.currentTimeMillis();
		LimitDeterministicAutomaton<EquivalenceClass, GeneralizedBreakpointFreeState, GeneralizedBuchiAcceptance, FGObligations> ldba =
				TRANSLATOR.apply(ltlFormula);
		logger.println("The automaton has " + ldba.size() + " states.");
		LDBAWrapper da = new LDBAWrapper(ldba);
		time = System.currentTimeMillis() - time;
		logger.println("Time for translation: " + (double) time / 1000.0 + " seconds.");

		return da;
	}

	public static int getInitialState()
	{
		return 0;
	}

	public int getAcceptanceSize()
	{
		return ldba.getAcceptingComponent().getAcceptance().size;
	}

	/**
	 * Get the destination of the edge from state i with label lab
	 *
	 * @param i     state
	 * @param label the label of the edge
	 * @return Returns -1 if no such edge is found.
	 */
	public int getEdgeDestByLabel(int i, NatBitSet label)
	{
		return get(getEdgeByLabel(i, label));
	}

	/**
	 * Get the destination of the edge from state i with label lab
	 *
	 * @param i     state
	 * @param label the label of the edge
	 * @return Returns null if no such edge is found.
	 */
	public Edge<?> getEdgeByLabel(int i, NatBitSet label)
	{
		Object state = get(i);

		if (state == null) {
			return null;
		}

		if (state instanceof EquivalenceClass) {
			return initialComponent.getEdge((EquivalenceClass) state, NatBitSets.toBitSet(label));
		}

		return acceptingComponent.getEdge((GeneralizedBreakpointFreeState) state, NatBitSets.toBitSet(label));
	}

	public int[] getJumpTargets(int i)
	{
		return jumpCache.computeIfAbsent(i, k -> {
			Object state = get(k);

			if (state instanceof EquivalenceClass) {
				Set<GeneralizedBreakpointFreeState> jumpTargets = ldba.getEpsilonJumps((EquivalenceClass) state);

				int[] targets = new int[jumpTargets.size()];
				int c = 0;

				for (GeneralizedBreakpointFreeState target : jumpTargets) {
					targets[c] = get(target);
					c++;
				}

				return targets;
			}

			return IntArrays.EMPTY_ARRAY;
		});
	}

	public int getInitialComponentSize()
	{
		return smallestAcceptingState;
	}

	public int size()
	{
		return ldba.size();
	}

	@Override
	public String toString()
	{
		return ldba.toString();
	}

	private Object get(int i)
	{
		return map.get(i);
	}

	public int get(Edge<?> edge)
	{
		return edge == null ? -1 : get(edge.getSuccessor());
	}

	public int get(Object state)
	{
		if (state == null) {
			return -1;
		}

		int i = numbering.getInt(state);

		if (i == -2) {
			i = numbering.size();
			map.add(state);
			numbering.put(state, i);
		}

		return i;
	}
}
