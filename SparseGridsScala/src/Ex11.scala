package demo

/*
 * This file contains an example solution to the Scala exercises for BPRD, by
 * David Christiansen, drc@itu.dk.
 */

/*
 * Exercises 11.1 through 11.4
 */
object Ex11 extends App {

  sealed abstract class Expr
  case class CstI(value: Int) extends Expr
  case class Var(name: String) extends Expr // 11.1 (ii)
  case class Prim(op: String, e1: Expr, e2: Expr) extends Expr

  /*
   * The compiler will give an incomplete pattern match warning below due to the
   * addition of Var in 11.1 (ii).
   */
  def eval0(expr: Expr): Int =
    expr match {
      case CstI(i) => i
      case Prim(op, e1, e2) =>
        val v1 = eval0(e1)
        val v2 = eval0(e2)
        op match {
          case "+" => v1 + v2
          case "-" => v1 - v2
          case "*" => v1 * v2
          case "/" => v1 / v2      // 11.1 (i)
          case "<" if v1 < v2 => 1 // 11.1 (i)
          case "<"            => 0 // 11.1 (i)
          case ">" if v1 > v2 => 1 // 11.1 (i)
          case ">"            => 0 // 11.1 (i)
          case "&" if v1 != 0 && v2 != 0 => 1 // 11.1 (i)
          case "&"                       => 0 // 11.1 (i)
        }
    }

  def eval1(expr: Expr, env: Map[String, Int]): Int =
    expr match {
      case CstI(i) => i
      case Var(x) => env(x)  // throws exception if x not in env
      case Prim(op, e1, e2) =>
        val v1 = eval1(e1, env)
        val v2 = eval1(e2, env)
        op match {
          case "+" => v1 + v2
          case "-" => v1 - v2
          case "*" => v1 * v2
          case "/" => v1 / v2      // 11.1 (i)
          case "<" if v1 < v2 => 1 // 11.1 (i)
          case "<"            => 0 // 11.1 (i)
          case ">" if v1 > v2 => 1 // 11.1 (i)
          case ">"            => 0 // 11.1 (i)
          case "&" if v1 != 0 && v2 != 0 => 1 // 11.1 (i)
          case "&"                       => 0 // 11.1 (i)
        }
    }

  val ex1 = Prim("*", CstI(7), Prim("+", Prim("-", CstI(10), CstI(6)), CstI(5)))
  val ex2 = Prim("*", CstI(7), Prim("+", Prim("<", CstI(10), CstI(6)), CstI(5)))
  val exVar = Prim("+", Var("x"), Var("y")) // 11.1 (ii)
  val testEnv = Map("x" -> 27, "y" -> 42)   // 11.1 (ii)

  println("Example 1: " + ex1.toString)
  println("         = " + eval0(ex1).toString)                // 11.1 (i)
  println("Example 2: " + ex2.toString)
  println("         = " + eval0(ex2).toString)                // 11.1 (i)
  println("Example with variables: " + exVar.toString + " in " + testEnv.toString) // 11.1 (ii)
  println("                      = " + eval1(exVar, testEnv).toString)             // 11.1 (ii)

  /*
   * Solution to 11.2.
   *
   * This implementation is one of many correct answers. I like it for the
   * following reasons:
   *  - it is easy to see the correspondence with the table in the exercise
   *  - subtrees are simplified, leading to more total simplification
   *  - the logic is primarily on the left of the pattern match rather than
   *    being hidden in if-expressions on the right of the pattern match. This
   *    is straightforward to read
   *  - nested pattern matches are avoided
   *  - the solution is concise
   *  Very similar code can be written in F#.
   */
  def simplify(expr: Expr): Expr =
    expr match {
      case Prim(op, e1, e2) =>
        (op, simplify(e1), simplify(e2)) match {
          case ("+", CstI(0), e      ) => e
          case ("+", e,       CstI(0)) => e
          case ("-", e,       CstI(0)) => e
          case ("*", CstI(1), e      ) => e
          case ("*", e,       CstI(1)) => e
          case ("*", CstI(0), _      ) => CstI(0)
          case ("*", _,       CstI(0)) => CstI(0)
          case ("-", e1s,     e2s    ) if e1s == e2s => CstI(0)
          case (op,  e1s,     e2s    ) => Prim(op, e1s, e2s)
        }
      case e => e // variables and constants are already simple
    }

  /*
   * Solution to 11.3.
   *
   * This implementation is ugly, as many of you noticed.  This was partly the
   * point, to show how much the for-expression can improve it. However, it also
   * has advantages over some of the alternatives.
   *
   * If we used the isEmpty and get methods, then we would need to always make
   * sure that we checked for isEmpty before calling get.  In other words, we'd
   * be just as badly off as if we just used null. The point of an option type
   * is that the type system should force us to do an empty check - isEmpty and
   * get are included as a sort of ugly backdoor for special circumstances.
   *
   * The pattern match ensures that the check for emptiness has taken place.
   */
  def eval2(expr: Expr, env: Map[String, Int]): Option[Int] =
    expr match {
      case CstI(i) => Some(i)
      case Var(x) => env get x // returns Some if x is in env, or None if not
      case Prim(op, e1, e2) =>
        eval2(e1, env) match {
          case None => None
          case Some(v1) =>
            eval2(e2, env) match {
              case None => None
              case Some(v2) =>
                op match {
                  case "+"            => Some(v1 + v2)
                  case "-"            => Some(v1 - v2)
                  case "*"            => Some(v1 * v2)
                  case "/" if v2 == 0 => None
                  case "/"            => Some(v1 / v2)
                  case "<" if v1 < v2 => Some(1)
                  case "<"            => Some(0)
                  case ">" if v1 > v2 => Some(1)
                  case ">"            => Some(0)
                  case "&" if v1 != 0 && v2 != 0 => Some(1)
                  case "&"                       => Some(0)
                }
            }
        }

    }

  /*
   * This is the solution to 11.4 that Hint 2 should lead you to.
   *
   * The reason that curly braces are used with for is to have Scala insert
   * the semicolons for us. Parentheses require explicit semicolons between
   * bindings.
   */
  def eval3(expr: Expr, env: Map[String, Int]): Option[Int] =
    expr match {
      case CstI(i) => Some(i)
      case Var(x) => env get x // returns Some if x is in env, or None if not
      case Prim(op, e1, e2) =>
        for {
          v1 <- eval3(e1, env)
          v2 <- eval3(e2, env)
          res <- op match {
                   case "+"            => Some(v1 + v2)
                   case "-"            => Some(v1 - v2)
                   case "*"            => Some(v1 * v2)
                   case "/" if v2 == 0 => None
                   case "/"            => Some(v1 / v2)
                   case "<" if v1 < v2 => Some(1)
                   case "<"            => Some(0)
                   case ">" if v1 > v2 => Some(1)
                   case ">"            => Some(0)
                   case "&" if v1 != 0 && v2 != 0 => Some(1)
                   case "&"                       => Some(0)
                 }
        } yield res
    }

  /*
   * This is an alternative way of solving 11.4 that you may think is prettier.
   *
   * For-expressions are allowed to contain filters in addition to bindings. In
   * the case of Option, a filter causes the expression to return None if the
   * condition is not fulfilled, or if it would be None anyway. If the condition
   * is fulfilled, then Some is left alone.
   *
   * This solution has a few advantages:
   *  - it is not necessary to write Some on the right-hand side of the pattern
   *    match, because a for ... yield expression causes the result to be
   *    automatically "packaged up" in the right constructor
   *  - it allows the evaluator to be more easily changed to use a different
   *    result type, in the manner shown in BPRD Lecture 12 on monads
   *  - all failure conditions are collected in the first part of the for, which
   *    may make the code more readable
   *
   *  There's also a disadvantage: the failure case on / is moved further from
   *  the place where the correct evaluation of / occurs. This may make the code
   *  less readable.
   */
  def eval3alt(expr: Expr, env: Map[String, Int]): Option[Int] =
    expr match {
      case CstI(i) => Some(i)
      case Var(x) => env get x // returns Some if x is in env, or None if not
      case Prim(op, e1, e2) =>
        for {
          v1 <- eval3alt(e1, env)
          v2 <- eval3alt(e2, env)
          if op != "/" || v2 != 0
        } yield (op match {
                   case "+"            => v1 + v2
                   case "-"            => v1 - v2
                   case "*"            => v1 * v2
                   case "/"            => v1 / v2
                   case "<" if v1 < v2 => 1
                   case "<"            => 0
                   case ">" if v1 > v2 => 1
                   case ">"            => 0
                   case "&" if v1 != 0 && v2 != 0 => 1
                   case "&"                       => 0
                 })
    }

}

/* Exercise 11.5 */
object Ex115 {
  sealed abstract class TypedExpr[T]
  //case class CstI(n: Int) extends TypedExpr[Int]
  //case class CstB(b: Boolean) extends TypedExpr[Boolean] // 11.5 (ii)
  case class Const[T](x: T) extends TypedExpr[T] // 11.5 (iii)

  case class Plus(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Int]
  case class Minus(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Int] // 11.5 (ii)
  case class Times(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Int]

  case class LessThanEq(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean]
  case class GreaterThanEq(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean] // 11.5 (ii)
  case class LessThan(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean] // 11.5 (ii)
  case class GreaterThan(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean] // 11.5 (ii)
  case class Equal(e1: TypedExpr[Int], e2: TypedExpr[Int]) extends TypedExpr[Boolean] // 11.5 (ii)

  case class IfThenElse[T](cond: TypedExpr[Boolean], e1: TypedExpr[T], e2: TypedExpr[T]) extends TypedExpr[T]

  /*
   * Many of you had a difficult time on 11.5 (iii).  In particular, many of you
   * attempted to inspect the runtime type of the contents of Const, which as
   * you know is difficult on the JVM in the presence of type erasure.
   *
   * The simple solution is to write that case of the pattern match
   * parametrically, which is to say not inspecting the contents of Const but
   * simply returning them.  This allows the compiler to check the types.
   *
   * This parametricity is exactly the same condition that allows F# to create
   * a generic function, and when you can use it, you should, even in languages
   * like C# where it's easier to inspect runtime types and values, because
   * it can lead to much cleaner code.
   */
  def eval[T](e: TypedExpr[T]): T = e match {
    case Const(x) => x // 11.5 (iii)
    //case CstI(n) => n
    //case CstB(b) => b // 11.5 (ii)
    case Plus(e1, e2) => eval(e1) + eval(e2)
    case Minus(e1, e2) => eval(e1) - eval(e2) // 11.5 (ii)
    case Times(e1, e2) => eval(e1) * eval(e2)
    case LessThanEq(e1, e2) => eval(e1) <= eval(e2)
    case GreaterThanEq(e1, e2) => eval(e1) >= eval(e2) // 11.5 (ii)
    case LessThan(e1, e2) => eval(e1) < eval(e2) // 11.5 (ii)
    case GreaterThan(e1, e2) => eval(e1) > eval(e2) // 11.5 (ii)
    case Equal(e1, e2) => eval(e1) == eval(e2) // 11.5 (ii)
    case IfThenElse(cond, e1, e2) => if (eval(cond)) eval(e1) else eval(e2)
  }

  // Begin 11.5 (i)
  val ex1 = Plus(Const(2), Times(Const(13), Const(4)))
  val ex2 = IfThenElse(LessThanEq(Const(23), Const(17)),
                       Const(7),
                       Plus(Const(3), Const(21)))
  // Fails because condition must be Boolean
  // val ex3 = IfThenElse(Const(0), Const(1), Const(2))
  // Fails because both branches of IfThenElse must have same type
  // val ex4 = IfThenElse(LessThanEq(Const(12), Const(3)), Const(42), LessThanEq(Const(3), Const(7)))
  val ex5 = IfThenElse(LessThanEq(Const(12), Const(3)),
                       LessThanEq(Const(42), Const(100)),
                       LessThanEq(Const(3), Const(7)))
  // End 11.5 (i)

}

