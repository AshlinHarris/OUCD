# OUCD
For orthogonal king graph tilings, produce an Optimal Unique Complete Dominating tile.

Consider a chess board with infinite length and width.
Each square has 8 neighbors: 4 on the sides, and 4 on the diagonals.
If we place a game piece on any square, we say that it "dominates" the 8 neighboring squares
(In our version of the rules, we say it doesn't dominates the square it covers).
If we assign a name to each piece, we can label every square by the set of pieces that dominate it.
We say the board is "completely" dominated if each square is next to at least one piece.
We say the board is "uniquely" dominated if no two squares have the same labels.
What fraction of squares must contain a piece in order for the entire board to be completely and uniquely dominated?

Solutions exist that cover 1/4 of all squares.
The simplest is the following 4 x 4 tile:
```
┼───┼───┼───┼───┼
│ @ │   │   │   │
┼───┼───┼───┼───┼
│   │ @ │   │   │
┼───┼───┼───┼───┼
│   │   │   │ @ │
┼───┼───┼───┼───┼
│   │   │ @ │   │
┼───┼───┼───┼───┼
```
To my knowledge, no one has proved that this is a lower bound.

This program generates finite rectangular boards to be tiled as a grid.
For any given size, it finds an optimal tile (that is, it uses the fewest number of pieces).
