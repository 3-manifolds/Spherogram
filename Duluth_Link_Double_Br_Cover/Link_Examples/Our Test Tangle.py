TestTangle = BraidTangle([2, 2, 1, 1, 2, 2, 1, 1, 2]) + IdentityBraid(1)
TestTangle.annular_closure().exterior().identify()
# [v3383(0,0)(0,0)(0,0), 8^3_1(0,0)(0,0)(0,0), L8a18(0,0)(0,0)(0,0)]
