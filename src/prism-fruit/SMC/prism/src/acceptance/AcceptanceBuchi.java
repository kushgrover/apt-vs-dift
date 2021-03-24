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

package acceptance;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import jdd.JDDVars;
import prism.PrismException;
import prism.PrismNotSupportedException;

import java.io.PrintStream;

/**
 * A Büchi acceptance condition (based on NatBitSet state sets).
 * The acceptance is defined via a set of "accepting" states
 * (sometimes also called final states) and is accepting if
 *  "infinitely often an accepting state is visited"
 */
public class AcceptanceBuchi implements AcceptanceOmega
{
	/** The set of goal states */
	private NatBitSet acceptingStates = NatBitSets.set();

	/** Constructor (no accepting states) */
	public AcceptanceBuchi()
	{
	}

	/** Constructor (set accepting states) */
	public AcceptanceBuchi(NatBitSet acceptingStates)
	{
		this.acceptingStates = acceptingStates;
	}

	/** Get the accepting state set */
	public NatBitSet getAcceptingStates()
	{
		return acceptingStates;
	}

	/** Set the accepting state set */
	public void setAcceptingStates(NatBitSet acceptingStates)
	{
		this.acceptingStates = acceptingStates;
	}

	/** Make a copy of the acceptance condition. */
	@Override public AcceptanceBuchi clone()
	{
		return new AcceptanceBuchi(acceptingStates.clone());
	}

	@Override
	public boolean isBSCCAccepting(NatBitSet bscc_states)
	{
		return bscc_states.intersects(acceptingStates);
	}

	/**
	 * Get the Rabin acceptance condition that is the equivalent of this Buchi condition.
	 */
	public AcceptanceRabin toRabin(int numStates)
	{
		AcceptanceRabin rabin = new AcceptanceRabin();
		rabin.add(new AcceptanceRabin.RabinPair(NatBitSets.set(), acceptingStates.clone()));
		return rabin;
	}

	/**
	 * Get the Streett acceptance condition that is the equivalent of this Buchi condition.
	 */
	public AcceptanceStreett toStreett(int numStates)
	{
		AcceptanceStreett streett = new AcceptanceStreett();
		NatBitSet allStates = NatBitSets.fullSet(numStates);
		streett.add(new AcceptanceStreett.StreettPair(allStates, acceptingStates.clone()));
		return streett;
	}

	/**
	 * Get a Rabin acceptance condition that is the complement of this condition, i.e.,
	 * any word that is accepted by this condition is rejected by the returned Rabin condition.
	 *
	 * @param numStates the number of states in the underlying model / automaton (needed for complementing NatBitSets)
	 * @return the complement Rabin acceptance condition
	 */
	public AcceptanceRabin complementToRabin(int numStates)
	{
		AcceptanceRabin rabin = new AcceptanceRabin();
		NatBitSet allStates = NatBitSets.fullSet(numStates);
		rabin.add(new AcceptanceRabin.RabinPair(acceptingStates.clone(), allStates));
		return rabin;
	}

	/**
	 * Get a Streett acceptance condition that is the complement of this condition, i.e.,
	 * any word that is accepted by this condition is rejected by the returned Streett condition.
	 * <br>
	 * Relies on the fact that once the goal states have been reached, all subsequent states
	 * are goal states.
	 *
	 * @param numStates the number of states in the underlying model / automaton (needed for complementing NatBitSets)
	 * @return the complement Streett acceptance condition
	 */
	public AcceptanceStreett complementToStreett(int numStates)
	{
		AcceptanceStreett streett = new AcceptanceStreett();
		streett.add(new AcceptanceStreett.StreettPair(acceptingStates.clone(), NatBitSets.set()));
		return streett;
	}

	/** Complement this acceptance condition, return as AcceptanceGeneric. */
	public AcceptanceGeneric complementToGeneric()
	{
		return toAcceptanceGeneric().complementToGeneric();
	}

	@Override
	public AcceptanceOmega complement(int numStates, AcceptanceType... allowedAcceptance) throws PrismException
	{
		if (AcceptanceType.contains(allowedAcceptance, AcceptanceType.RABIN)) {
			return complementToRabin(numStates);
		} else if (AcceptanceType.contains(allowedAcceptance, AcceptanceType.STREETT)) {
			return complementToStreett(numStates);
		} else if (AcceptanceType.contains(allowedAcceptance, AcceptanceType.GENERIC)) {
			return complementToGeneric();
		}
		throw new PrismNotSupportedException("Can not complement " + getType() + " acceptance to a supported acceptance type");
	}

	@Override
	public void lift(LiftNatBitSet lifter)
	{
		acceptingStates=lifter.lift(acceptingStates);
	}

	@Override
	public AcceptanceBuchiDD toAcceptanceDD(JDDVars ddRowVars)
	{
		return new AcceptanceBuchiDD(this, ddRowVars);
	}

	@Override
	public AcceptanceGeneric toAcceptanceGeneric()
	{
		return new AcceptanceGeneric(AcceptanceGeneric.ElementType.INF, acceptingStates.clone());
	}

	@Override
	public String getSignatureForState(int i)
	{
		return acceptingStates.contains(i) ? "!" : " ";
	}

	@Override
	public String getSignatureForStateHOA(int stateIndex)
	{
		if (acceptingStates.contains(stateIndex)) {
			return "{0}";
		} else {
			return "";
		}
	}

	/** Returns a textual representation of this acceptance condition. */
	@Override
	public String toString()
	{
		return acceptingStates.toString();
	}

	@Override
	public String getSizeStatistics()
	{
		return acceptingStates.size() + " accepting states";
	}

	@Override
	public AcceptanceType getType()
	{
		return AcceptanceType.BUCHI;
	}

	@Override
	@Deprecated
	public String getTypeAbbreviated()
	{
		return getType().getNameAbbreviated();
	}

	@Override
	@Deprecated
	public String getTypeName()
	{
		return getType().getName();
	}

	@Override
	public void outputHOAHeader(PrintStream out)
	{
		out.println("acc-name: Buchi");
		out.println("Acceptance: 1 Inf(0)");
	}
}
