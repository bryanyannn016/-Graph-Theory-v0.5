/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graphtheory;

import java.awt.Color;
import java.awt.Graphics;
import javax.swing.JOptionPane;
/**
 *
 * @author mk
 */
public class Edge {

    public Vertex vertex1;
    public Vertex vertex2;
    public boolean wasFocused;
    public boolean wasClicked;
    public Integer weight; // New weight attribute

    public Edge(Vertex v1, Vertex v2) {
        this.vertex1 = v1;
        this.vertex2 = v2;
        this.weight = promptForWeight(); // Prompt user for weight
    }

    // Prompt user for edge weight
    private Integer promptForWeight() {
        String[] options = {"Enter Weight", "No Weight", "Cancel"};
        int choice = JOptionPane.showOptionDialog(
            null,
            "Choose an option for edge weight:",
            "Edge Weight",
            JOptionPane.DEFAULT_OPTION,
            JOptionPane.PLAIN_MESSAGE,
            null,
            options,
            options[0]
        );

        if (choice == 0) { // Enter Weight
            while (true) {
                String input = JOptionPane.showInputDialog(null, "Enter edge weight:", "Edge Weight", JOptionPane.PLAIN_MESSAGE);
                if (input == null) { // User canceled
                    return 1; // Default weight
                }
                try {
                    return Integer.parseInt(input);
                } catch (NumberFormatException e) {
                    JOptionPane.showMessageDialog(null, "Invalid input! Please enter a number.");
                }
            }
        } else if (choice == 1) { // No Weight
            return null; // No weight assigned
        }

        return 1; // Default weight if canceled
    }

    public void draw(Graphics g) {
        if (wasClicked) {
            g.setColor(Color.red);
        } else if (wasFocused) {
            g.setColor(Color.blue);
        } else {
            g.setColor(Color.black);
        }

        int x1 = vertex1.location.x;
        int y1 = vertex1.location.y;
        int x2 = vertex2.location.x;
        int y2 = vertex2.location.y;

        g.drawLine(x1, y1, x2, y2);

        // Draw weight in the middle of the edge if it has one
        if (weight != null) {
            int midX = (x1 + x2) / 2;
            int midY = (y1 + y2) / 2;
            g.setColor(Color.BLACK);
            g.drawString(String.valueOf(weight), midX, midY);
        }
    }


    public boolean hasIntersection(int x, int y) {
        int x1, x2, y1, y2;
        x1 = vertex1.location.x;
        x2 = vertex2.location.x;
        y1 = vertex1.location.y;
        y2 = vertex2.location.y;
        float slope = 0;
        if (x2 != x1) {
            slope = (y2 - y1) / (x2 - x1);
        }

        float b = Math.abs(x1 * slope - y1);

        if (y + b <= Math.round(slope * x) + 10 && y + b >= Math.round(slope * x) - 10) {
            if (x1 > x2 && y1 > y2) {
                if (x <= x1 && x >= x2 && y <= y1 && y >= y2) {
                    return true;
                }
            } else if (x1 < x2 && y1 > y2) {
                if (x <= x2 && x >= x1 && y <= y1 && y >= y2) {
                    return true;
                }
            } else if (x1 < x2 && y1 < y2) {
                if (x <= x2 && x >= x1 && y <= y2 && y >= y1) {
                    return true;
                }
            } else if (x <= x1 && x >= x2 && y <= y2 && y >= y1) {
                return true;
            }
        }
        return false;

    }
}
