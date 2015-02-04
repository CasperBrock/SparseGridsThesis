//BPRD Assignment 11 by Casper Brock (casb@itu.dk)
//Exercise 11.5

sealed abstract class TypedExpr[T]
//case class CstI(n: Int) extends TypedExpr[Int]
case class Plus(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Int]
case class Times(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Int]
case class LessThanEq(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean]
case class IfThenElse[T](cond: TypedExpr[Boolean], e1: TypedExpr[T], e2: TypedExpr[T]) extends TypedExpr[T]

//Exercise 11.5 (ii)
//case class CstB(n: Boolean) extends TypedExpr[Boolean]
case class Sub(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Int]
case class LessThan(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean]
case class GreaterThan(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean]
case class GreaterThanEq(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean]
case class Eq(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean]
case class And(e1: TypedExpr[Boolean], e2: TypedExpr[Boolean]) extends TypedExpr[Boolean]
case class Or(e1: TypedExpr[Boolean], e2: TypedExpr[Boolean]) extends TypedExpr[Boolean]

//Exercise 11.5 (iii)
case class Const[T] (n: T) extends TypedExpr[T]

object exercise112 {
def main(args: Array[String]) = { 

  //Exercise 11.5 (i)
  val v1 = Plus(CstI(2), Times(CstI(13), CstI(4)))
  val v2 = IfThenElse(LessThanEq(CstI(23), CstI(17)), CstI(7), Plus(CstI(3), CstI(21)))
  //val v3 = IfThenElse(CstI(0), CstI(1), CstI(2))  //Doesn't type check, e1 must have type Boolean
  //val v4 = IfThenElse(LessThanEq(CstI(12), CstI(3)), CstI(42), LessThanEq(CstI(3), CstI(7)))  //Doesn't type check (e2 and e3 must have same type)
  val v5 = IfThenElse(LessThanEq(CstI(12), CstI(3)), LessThanEq(CstI(42), CstI(1000)), LessThanEq(CstI(3), CstI(7)))
}

def eval[T](e: TypedExpr[T]): T = 
  e match {
    //case CstI(n) => n
    case Plus(e1, e2) => eval(e1) + eval(e2)
    case Times(e1, e2) => eval(e1) * eval(e2)
    case LessThanEq(e1, e2) => eval(e1) <= eval(e2)
    case IfThenElse(cond, e1, e2) => if (eval(cond)) eval(e1) else eval(e2)
    
    //Exercise 11.5 (ii)
    //case CstB(n) => n
    case Sub(e1, e2) => eval(e1) - eval(e2)
    case LessThan(e1, e2) => eval(e1) < eval(e2)
    case GreaterThanEq(e1, e2) => eval(e1) >= eval(e2)
    case GreaterThan(e1, e2) => eval(e1) > eval(e2)
    case Eq(e1, e2) => eval(e1) == eval(e2)
    case And(e1, e2) => eval(e1) && eval(e2)
    case Or(e1, e2) => eval(e1) || eval(e2)

    //Exercise 11.5 (iii)
    case Const(n) => n
  }
}
