object gg {

	def fact(N: Int): Int = 
    N match {
      case N if N > 0 => fact(N - 1)
      case 1 => N
	}
}