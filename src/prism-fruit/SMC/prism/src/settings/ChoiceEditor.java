//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Andrew Hinton <ug60axh@cs.bham.ac.uk> (University of Birmingham)
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

package settings;

import javax.swing.JComboBox;
import javax.swing.JTable;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.util.List;

public class ChoiceEditor implements SettingEditor, ActionListener, FocusListener
{
	private Font font = new Font("monospaced", Font.PLAIN, 12);
	private String[] choices;
	private JComboBox<String> combo;
	private JTable lastTable;
	private boolean modified = false;

	/** Creates a new instance of ChoiceEditor */
	public ChoiceEditor(String[] choices)
	{
		combo = new JComboBox<>(choices);
		combo.addActionListener(this);
		combo.addFocusListener(this);
		this.choices = choices;
	}

	@Override public Object getEditorValue() throws SettingException
	{
		if (modified) {
			modified = false;
			if (combo.getSelectedItem() != null)
				return combo.getSelectedItem().toString();
			else
				return NOT_CHANGED_VALUE;
		} else
			return NOT_CHANGED_VALUE;
	}

	@Override public Component getTableCellEditorComponent(JTable table, Setting owner, Object value, boolean isSelected, int row, int column)
	{

		if (isSelected) {
			combo.setForeground(table.getSelectionForeground());
			combo.setBackground(table.getSelectionBackground());
		} else {
			combo.setForeground(table.getForeground());
			combo.setBackground(table.getBackground());
		}

		combo.setBorder(new EmptyBorder(0, 0, 0, 0));//UIManager.getBorder("Table.focusCellHighlightBorder") );
		combo.setFont(font);
		combo.setFocusable(false);

		if (value instanceof String) {
			int index = -1;
			for (int i = 0; i < choices.length; i++) {
				if (choices[i].equals(value.toString())) {
					index = i;
					break;
				}
			}
			combo.setSelectedIndex(index);
		} else if (value instanceof List) {
			List<?> values = (List<?>) value;
			if (!values.isEmpty()) {
				//if we have multiple properties selected.
				String last = null;
				boolean allSame = true;
				for (Object val : values) {
					if (val instanceof String) {
						String str = (String) val;
						if (last != null) {
							if (!str.equals(last)) {
								allSame = false;
								break;
							}
							last = str;
						} else {
							last = str;
						}
					}
				}
				if (allSame) {

					int index = -1;
					for (int i = 0; i < choices.length; i++) {
						if (choices[i].equals(last)) {
							index = i;
							break;
						}
					}
					combo.setSelectedIndex(index);
				} else {
					combo.setSelectedIndex(-1);

				}

			}
		}

		lastTable = table;
		return combo;
	}

	@Override public void stopEditing()
	{
	}

	@Override public void actionPerformed(ActionEvent e)
	{
		modified = true;
		if (lastTable != null)
			lastTable.editingStopped(new ChangeEvent(this));
	}

	@Override public void focusGained(FocusEvent e)
	{
	}

	@Override public void focusLost(FocusEvent e)
	{
		if (lastTable.getCellEditor() != null)
			lastTable.removeEditor();
	}

}
