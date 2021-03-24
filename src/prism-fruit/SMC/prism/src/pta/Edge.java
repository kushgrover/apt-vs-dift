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

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import prism.PrismException;

import java.io.IOException;
import java.io.Writer;

public class Edge
{
	// Parent transition;
	private Transition parent;
	// Probability
	private double prob;
	// Destination location
	private int dest;
	// Resets
	private final Int2IntMap resets;

	/**
	 * Create an empty edge.
	 * @param parent Parent transition
	 * @param prob Probability
	 * @param dest destination location
	 */
	public Edge(Transition parent, double prob, int dest)
	{
		this.parent = parent;
		this.prob = prob;
		this.dest = dest;
		resets = new Int2IntOpenHashMap();
	}

	/**
	 * Copy constructor.
	 */
	public Edge(Edge edge)
	{
		this(edge.parent, edge.prob, edge.dest);
		for (Int2IntMap.Entry e : edge.resets.int2IntEntrySet()) {
			addReset(e.getIntKey(), e.getIntValue());
		}
	}

	public void setProb(double prob)
	{
		this.prob = prob;
	}

	public void setDestination(int dest)
	{
		this.dest = dest;
	}

	public void setParent(Transition parent)
	{
		this.parent = parent;
	}

	public void addReset(int clock)
	{
		addReset(clock, 0);
	}

	public void addReset(int clock, int val)
	{
		resets.put(clock, val);
	}

	public Transition getParent()
	{
		return parent;
	}

	public int getDestination()
	{
		return dest;
	}

	public double getProbability()
	{
		return prob;
	}

	public Iterable<Int2IntMap.Entry> getResets()
	{
		return resets.int2IntEntrySet();
	}

	/**
	 * Perform some basic syntactic checks.
	 */
	public void check() throws PrismException
	{
	}

	public String toString()
	{
		boolean first = true;
		StringBuilder s = new StringBuilder();
		s.append(prob).append(" : {");
		for (Int2IntMap.Entry e : resets.int2IntEntrySet()) {
			if (first)
				first = false;
			else
				s.append(",");
			s.append(parent.getParent().getClockName(e.getIntKey()));
			s.append("=").append(e.getIntValue());
		}
		s.append("}").append(parent.getParent().getLocationName(dest));
		return s.toString();
	}

	public void writeToDesFile(Writer out, String actionPrefix, Iterable<Constraint> guard) throws PrismException,
			IOException
	{
		boolean first;
		PTA pta = parent.getParent();
		out.write("\t" + actionPrefix + "tran  ");
		out.write(pta.getLocationName(dest).toString().replace(':', '_') + "; ");
		out.write(Constraint.toStringList(pta, guard));
		out.write("; ");
		first = true;
		for (Int2IntMap.Entry e : resets.int2IntEntrySet()) {
			if (first)
				first = false;
			else
				out.write(",");
			out.write(pta.getClockName(e.getIntKey()) + "=" + e.getIntValue());
		}
		if (first)
			out.write("null");
		out.write("; " + prob + "\n");
	}
}
